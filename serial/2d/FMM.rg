
--
-- Main driver of Fast Multipole Method simulation code
-- Authors : Brian C. Munguia
--

-- Import libraries
import "regent"
local c       = regentlib.c
local std     = terralib.includec("stdlib.h")
local sqrt    = regentlib.sqrt(double)
local cmath   = terralib.includec("math.h")
local cstring = terralib.includec("string.h")

-- Helper module
local FMMConfig = require("FMM_config")

-- 2D vector type
struct Vector2d {
  _1 : double;
  _2 : double;
}

-- Fieldspace 'Particle'
fspace Particle {
  q        : double;   -- Particle source strength
  pos      : Vector2d; -- Position vector
  vel      : Vector2d; -- Velocity vector
  force    : Vector2d; -- Force vector
  field    : Vector2d; -- Field strength
  boxes    : int2d[5]; -- Boxes containing particle
  id       : int32;    -- Particle ID
}

-- Fieldspace 'Box'
fspace Box {
  Parent     : int2d;      -- Parent box
  Children   : int2d[4];   -- Children boxes, up to 8 per box in 2d
  Ilist      : int2d[27]; -- Interaction list, up to 189 (6^2 - 3^2)
  neighb     : int2d[8];  -- Neighbors
  part       : int32[64];  -- Particles in box
  num_child  : int32;      -- Number of children
  num_ilist  : int32;      -- Number of interactions
  num_neighb : int32;      -- Number of neighbors
  num_part   : int32;      -- Number of particles in box
  ctr        : Vector2d;   -- Box center
  S          : double;     -- Box side length
  --phi        : Vector2d;   -- p-term multipole expansion due to particles in box
  --psi        : Vector2d;   -- p-term expansion due to interaction list
  --psi_tilde  : Vector2d;   -- p-term local expansion due to particles in parent's
                           -- interaction list
  a_0        : double;     -- Coefficient in multipole expansion
  a_k        : Vector2d[20]; -- Coefficient in multipole expansion
  b_l        : Vector2d[20]; -- Coefficient in local expansion
}

terra skip_header(f : &c.FILE)
  var x0 : int8[512], y0 : uint32, x1 : int8[512], y1 : int8[512]
  var x2 : int8[512], y2 : double
  c.fscanf(f, "%s %llu\n %s %s\n %s %lf", &x0, &y0, &x1, &y1, &x2, &y2)
end

terra get_force_const(model_type : &int8)
  var force_const : double
  if cstring.strcmp(model_type, "GRAVITY") == 0 then
    force_const = 6.6740831e-11
  elseif cstring.strcmp(model_type, "COULOMB") == 0 then
    force_const = -8.987551787e9
  else
    force_const = 1.0
  end
  return force_const
end

terra read_particles(f : &c.FILE, source : &double)
  var x : int8[512]
  return c.fscanf(f, "%s %lf %lf %lf %lf %lf %lf %lf\n", &x, &source[0],
                  &source[1], &source[2], &source[3], &source[4]) == 6
end

terra Vector2d:kfun()
  return 1.0/sqrt(self._1 * self._1 + self._2 * self._2)
end

terra Vector2d.metamethods.__add(vec1 : Vector2d, vec2 : Vector2d)
  return Vector2d { vec1._1 + vec2._1, vec1._2 + vec2._2 }
end

terra Vector2d.metamethods.__sub(vec1 : Vector2d, vec2 : Vector2d)
  return Vector2d { vec1._1 - vec2._1, vec1._2 - vec2._2 }
end

terra Vector2d.metamethods.__mul(c : double, v : Vector2d)
  return Vector2d { v._1 * c, v._2 * c }
end

terra well_separated(S : double, r : Vector2d)
  r._1 = cmath.round(sqrt(1.0/(S*S)) * r._1)
  r._2 = cmath.round(sqrt(1.0/(S*S)) * r._2)
  if r._1 > 1 or r._2 > 1 then
    return true
  else
    return false
  end
end

terra comp_mul(vec1 : Vector2d, vec2 : Vector2d)
  return Vector2d {vec1._1 * vec2._1 - vec1._2 * vec2._2,
                   vec1._1 * vec2._2 + vec1._2 * vec2._1}
end

terra comp_pow(vec : Vector2d, k : int32)
  var theta : double = cmath.atan2(vec._2, vec._1)
  var z : double = cmath.pow(sqrt(cmath.pow(vec._1,2)+cmath.pow(vec._1,2)),k)
  return Vector2d {z * cmath.cos(k*theta), z * cmath.sin(k*theta)}
end

terra comp_log(vec : Vector2d)
  return Vector2d {cmath.log(sqrt(cmath.pow(vec._1,2) + cmath.pow(vec._1,2))),
                   cmath.atan(vec._2/vec._1)}
end

terra factorial(n : uint32)
  var fact : uint32 = n
  if n == 0 then return 1
    else
    for k = 1, n+1 do
      fact = fact * k
    end
    return fact
  end
end

terra choose(n : uint32, k : uint32)
  return factorial(n)/(factorial(k) * factorial(n - k))
end

task initialize_particles(r_particles   : region(Particle),
                          filename      : int8[512])
where
  reads writes(r_particles)
do
  -- Set initial force and field to 0, and particle IDs
  var ipart : int32 = 0
  for particle in r_particles do
    particle.force = Vector2d {0.0, 0.0}
    particle.field = Vector2d {0.0, 0.0}
    particle.id = ipart
    ipart = ipart + 1
  end

  var f = c.fopen(filename, "rb")
  skip_header(f)

  -- Source data from input file (strength, position, velocity)
  var source : double[5]

  for particle in r_particles do
    regentlib.assert(read_particles(f, source), "Less data than it should be")
    particle.q = source[0]
    particle.pos = Vector2d {source[1], source[2]}
    particle.vel = Vector2d {source[3], source[4]}
  end

  c.fclose(f)

end

task create_tree(num_lvl     : uint32,
                 r_particles : region(Particle),
                 r_boxes     : region(ispace(int2d, box_per_dim_i2d), Box),
                 particle    : ptr(Particle,r_particles),
                 ctr         : Vector2d,
                 S           : double)
where
  reads writes(r_boxes, r_particles)
do
  var x : double = particle.pos._1
  var y : double = particle.pos._2

  particle.boxes[0] = {0,0}

  var icurr : int2d = {0, 0}
  var ibox : int32 = 0
  var jbox : int32 = 0
  var icell_old : int32 = 0
  var jcell_old : int32 = 0
  var i : int32 = 1
  var offset : int32 = -1
  for lvl = 0, num_lvl do
    -- Index of child cell containing point
    var offset : int32 = -1 + cmath.pow(2,lvl+1)
    var icell : uint32 = int(x > ctr._1)
    var jcell : uint32 = int(y > ctr._2)

    S = 0.5*S

    ctr._2 = ctr._2 + 0.5*S*(2.0*jcell-1.0)
    ctr._1 = ctr._1 + 0.5*S*(2.0*icell-1.0)

    -- Index of child node
    ibox = icell + 2*icell_old + offset
    jbox = jcell + 2*jcell_old + offset

    -- Check whether child node exists
    var child_exist : uint8 = 0
    for j = 0, 4 do
      if r_boxes[icurr].Children[j] == [int2d]{ibox,jbox} then
        child_exist = 1
        break
      end
    end

    if child_exist == 0 then
      for j = 0,4 do
        if r_boxes[icurr].Children[j] == [int2d]{-1,-1} then
        r_boxes[icurr].Children[j] = {ibox,jbox}
        r_boxes[icurr].num_child = r_boxes[icurr].num_child + 1
        r_boxes[r_boxes[icurr].Children[j]].Parent = icurr
        r_boxes[r_boxes[icurr].Children[j]].ctr = ctr
        r_boxes[r_boxes[icurr].Children[j]].S = S
        end
      end
    end

    r_boxes[{ibox,jbox}].num_part = r_boxes[{ibox,jbox}].num_part + 1
    particle.boxes[lvl+1] = {ibox,jbox}
    r_boxes[{ibox,jbox}].part[r_boxes[{ibox,jbox}].num_part-1] =
      unsafe_cast(ptr(Particle,r_particles),particle.id)

    -- Update current box
    icurr = {ibox,jbox}
    i = i+1
    icell_old = icell
    jcell_old = jcell

    var box = r_boxes[{ibox,jbox}]

  end

end

task create_neighb(r_boxes : region(ispace(int2d, box_per_dim_i2d), Box),
                   num_lvl : uint32)
where
  reads writes(r_boxes.num_neighb, r_boxes.neighb)
do
  for lvl = 1,num_lvl+1 do
  var box_per_dim_lvl : uint32 = cmath.pow(2,lvl)
  var offset : int32 = -1 + cmath.pow(2,lvl)
  for i = 0, box_per_dim_lvl do
  for j = 0, box_per_dim_lvl do
  -- Inner loop to find boxes within domain
    for ii = -1, 2 do
    for jj = -1, 2 do
      if (i+ii > -1) and (i+ii < box_per_dim_lvl) and (j+jj > -1) and
         (j+jj < box_per_dim_lvl) and (ii ~= 0 or jj ~= 0) then
           var num_neighb : int32 =
               r_boxes[{i+offset, j+offset}].num_neighb + 1
           r_boxes[{i+offset, j+offset}].num_neighb = num_neighb
           r_boxes[{i+offset, j+offset}].neighb[num_neighb-1] =
               {i+ii+offset, j+jj+offset}
      end
    end
    end
  end
  end
  end

end

task create_Ilist(r_boxes : region(ispace(int2d, box_per_dim_i2d), Box),
                  node_ind    : int2d)
where
  reads(r_boxes.num_neighb, r_boxes.neighb, r_boxes.num_ilist, r_boxes.Children,
        r_boxes.num_child, r_boxes.ctr, r_boxes.S, r_boxes.Parent, r_boxes.Ilist),
  writes(r_boxes.num_ilist, r_boxes.Ilist)
do
  var S : double = r_boxes[node_ind].S
  var ctr : Vector2d = r_boxes[node_ind].ctr

  if r_boxes[node_ind].Parent ~= [int2d]{-1,-1} then
    var num_pn : int32 = r_boxes[r_boxes[node_ind].Parent].num_neighb
    var pn : int2d[8] = r_boxes[r_boxes[node_ind].Parent].neighb
    -- Loop over parent's neighbors
    for j = 0, num_pn do
      -- Loop over children
      for i = 0, r_boxes[node_ind].num_child do
        var child : int2d = r_boxes[pn[j]].Children[i]
        var r : Vector2d = r_boxes[child].ctr - ctr
        if well_separated(S,r) then
          r_boxes[node_ind].num_ilist = r_boxes[node_ind].num_ilist + 1
          r_boxes[node_ind].Ilist[r_boxes[node_ind].num_ilist-1] = child
        end
      end
    end

  end

  for i = 0, 4 do
    if r_boxes[node_ind].Children[i] ~= [int2d]{-1,-1} then
      create_Ilist(r_boxes, r_boxes[node_ind].Children[i])
    end
  end

end

task init_Ilist_tree(r_particles : region(Particle),
                     r_boxes     : region(ispace(int2d, box_per_dim_i2d), Box),
                     ctr         : Vector2d,
                     S           : double,
                     num_lvl     : uint32,
                     p           : uint32)
where
  reads writes(r_particles, r_boxes)
do
  var i : uint32 = 0

  -- Initialize boxes
  r_boxes[{0,0}].ctr = ctr
  r_boxes[{0,0}].S = S
  var child_init : int2d[4]
  var a_k_init : Vector2d[20]
  for i = 0, 4 do
    child_init[i] = {-1, -1}
    a_k_init[i] = {0.0, 0.0}
  end
  for i = 4, p do
    a_k_init[i] = {0.0, 0.0}
  end
  fill(r_boxes.Children, child_init)
  fill(r_boxes.num_part, 0)
  fill(r_boxes.num_child, 0)
  fill(r_boxes.num_neighb, 0)
  fill(r_boxes.num_ilist, 0)
  fill(r_boxes.Parent, {-1, -1})
  fill(r_boxes.a_k, a_k_init)
  fill(r_boxes.b_l, a_k_init)
  -- Create FMM tree
  for particle in r_particles do
    create_tree(num_lvl, r_particles, r_boxes, particle, ctr, S)
  end

  -- Create neighbors list
  create_neighb(r_boxes, num_lvl)

  -- Create interactions list
  create_Ilist(r_boxes, int2d {0,0})

end

task P2M(r_particles : region(Particle),
         r_boxes     : region(ispace(int2d, box_per_dim_i2d), Box),
         particle    : ptr(Particle,r_particles),
         lvl         : uint32,
         p           : uint32)
where
  reads(r_particles.q, r_particles.pos, r_particles.boxes, r_boxes),
  writes(r_boxes.a_k)
do
  var ds : Vector2d = particle.pos - r_boxes[particle.boxes[lvl]].ctr

  if ds._1 ~= 0.0 or ds._2 ~= 0.0 then
    for k = 1, p do
      var a_k : Vector2d = -1/k * particle.q * comp_pow(ds,k)
      r_boxes[particle.boxes[lvl]].a_k[k] = r_boxes[particle.boxes[lvl]].a_k[k] + a_k
    end
    r_boxes[particle.boxes[lvl]].a_k[0] = r_boxes[particle.boxes[lvl]].a_k[0]
                                       + Vector2d {particle.q, 0.0}
  end

end

task outgoing_translation(r_boxes : region(ispace(int2d, box_per_dim_i2d), Box),
                          box_id  : int2d,
                          p       : uint32)
where
  reads(r_boxes), writes(r_boxes.a_k)
do
  for k = 0, p do
    if k == 0 then
      for j = 0, 4 do
        r_boxes[box_id].a_k[k] = r_boxes[box_id].a_k[k]
                                 + r_boxes[r_boxes[box_id].Children[j]].a_k[k]
      end
    else
      for j = 0, 4 do
        var d : Vector2d = r_boxes[box_id].ctr
                           - r_boxes[r_boxes[box_id].Children[j]].ctr
        var a_sig : Vector2d = r_boxes[r_boxes[box_id].Children[j]].a_k[0]
        r_boxes[box_id].a_k[k] = r_boxes[box_id].a_k[k]
                                 - 1/k * comp_mul(comp_pow(d,k),a_sig)
        for i = 1, k do
          var ch_pq : int32 = choose(k,i)
          a_sig = r_boxes[r_boxes[box_id].Children[j]].a_k[i]
          r_boxes[box_id].a_k[k] = r_boxes[box_id].a_k[k]
                                   + ch_pq * comp_mul(comp_pow(d, k - i), a_sig)
        end
      end
    end
  end
end

task incoming_translation(r_boxes : region(ispace(int2d, box_per_dim_i2d), Box),
                          box_id  : int2d,
                          p       : uint32)
where
  reads(r_boxes.num_ilist, r_boxes.Ilist, r_boxes.ctr, r_boxes.b_l, r_boxes.a_k),
  writes(r_boxes.b_l)
do
  if r_boxes[box_id].num_ilist > 0 then

  for l = 0, p do
    if l == 0 then
      for j = 0, r_boxes[box_id].num_ilist do
        var d : Vector2d =  r_boxes[r_boxes[box_id].Ilist[j]].ctr
                            - r_boxes[box_id].ctr
        r_boxes[box_id].b_l[l] = r_boxes[box_id].b_l[l]
                                 + comp_mul(r_boxes[r_boxes[box_id].Ilist[j]].a_k[0],
                                 comp_log(Vector2d {-d._1, -d._2}))
        for k = 1, p do
          r_boxes[box_id].b_l[l] = r_boxes[box_id].b_l[l]
                                   + cmath.pow(-1.0,k)
                                   * comp_mul(comp_pow(d, -k),
                                   r_boxes[r_boxes[box_id].Ilist[j]].a_k[k])
        end
      end

    else
      for j = 0, r_boxes[box_id].num_ilist do
        var d : Vector2d =  r_boxes[r_boxes[box_id].Ilist[j]].ctr
                            - r_boxes[box_id].ctr

        var b_new : Vector2d = -1/l * r_boxes[r_boxes[box_id].Ilist[j]].a_k[0]
        for k = 1, p do
          var ch_pq : int32 = choose(l+k-1,k-1)
          b_new = b_new + ch_pq * cmath.pow(-1,k)
                  * comp_mul(r_boxes[r_boxes[box_id].Ilist[j]].a_k[k], comp_pow(d,-k))
        end
        r_boxes[box_id].b_l[l] = r_boxes[box_id].b_l[l]
                                 + comp_mul(b_new, comp_pow(d,-l))

      end
    end
  end

  end
end

task L2L(r_boxes : region(ispace(int2d, box_per_dim_i2d), Box),
         box_id  : int2d,
         p       : uint32)
where
  reads(r_boxes), writes(r_boxes.b_l)
do

  for j = 0, r_boxes[box_id].num_child do
    var box_child = r_boxes[r_boxes[box_id].Children[j]]
    var d : Vector2d =  r_boxes[box_id].ctr - box_child.ctr
    box_child.b_l = r_boxes[box_id].b_l
    for l = 0, p do
      for k = p-j-1, p do
        var ch_pq : uint32 = choose(k,l)
        -- Horner scheme
        box_child.b_l[k] = box_child.b_l[k]
                           - comp_mul(d,r_boxes[box_id].b_l[k+1])
      end
    end
    r_boxes[r_boxes[box_id].Children[j]].b_l = box_child.b_l
  end

end

task M2M(r_particles : region(Particle),
         r_boxes     : region(ispace(int2d, box_per_dim_i2d), Box),
         num_lvl     : uint32,
         box_per_dim : uint32,
         p           : uint32)
where
  reads(r_particles, r_boxes), writes(r_boxes)
do
  -- Form multipole expansions of potential field due to particles in each box
  -- about the box center at the finest mesh level
  var lvl : uint32 = num_lvl
  for particle in r_particles do
    P2M(r_particles, r_boxes, particle, lvl, p)
  end

  -- Form multipole expansions about centers of all boxes at coarser mesh
  -- levels
  var box_min : uint32 = box_per_dim - cmath.pow(2,lvl)
  var box_max : uint32 = box_per_dim - 1
  for lvl = num_lvl-1, -1, -1 do
    box_min = box_min - cmath.pow(2,lvl)
    box_max = box_max - cmath.pow(2,lvl+1)
    for ibox = box_min, box_max+1 do
      for jbox = box_min, box_max+1 do
        outgoing_translation(r_boxes, {ibox, jbox}, p)
      end
    end
  end

  return 1

end

task M2L(r_boxes     : region(ispace(int2d, box_per_dim_i2d), Box),
         num_lvl     : uint32,
         box_per_dim : uint32,
         p           : uint32)
where
  reads(r_boxes), writes(r_boxes.b_l)
do
  var box_min : uint32 = 0
  var box_max : uint32 = 0
  for lvl = 1, num_lvl do
    box_min = box_min + cmath.pow(2,lvl-1)
    box_max = box_max + cmath.pow(2,lvl)
    c.printf("M2L : First inner loop\n")
    -- First inner loop converts multipole expansion of each box in interaction
    -- list to local expansion
    for ibox = box_min, box_max+1 do
      for jbox = box_min, box_max+1 do
        incoming_translation(r_boxes, {ibox, jbox}, p)
      end
    end

    c.printf("M2L : Second inner loop\n")
    -- Second inner loop forms local expansion by expanding from box about
    -- children box centers
    for ibox = box_min, box_max+1 do
      for jbox = box_min, box_max+1 do
        L2L(r_boxes, {ibox, jbox}, p)
      end
    end
  end

  c.printf("M2L : Finest mesh level\n")
  -- Compute interactions at finest mesh level
  box_min = box_min + cmath.pow(2,num_lvl-1)
  box_max = box_max + cmath.pow(2,num_lvl)
  for ibox = box_min, box_max+1 do
    for jbox = box_min, box_max+1 do
      incoming_translation(r_boxes, {ibox, jbox}, p)
    end
  end
  c.printf("M2L : complete\n")

  return 1

end

task eval_forces(r_particles : region(Particle),
                 r_boxes     : region(ispace(int2d, box_per_dim_i2d), Box),
                 num_lvl     : uint32,
                 box_per_dim : uint32,
                 p           : uint32)
where
  reads(r_particles, r_boxes),
  writes(r_particles.field, r_particles.force)
do

  for particle in r_particles do
    var box = r_boxes[particle.boxes[num_lvl]]
    var y = particle.pos
    var ctr = box.ctr
    particle.field = particle.field + box.b_l[0]--comp_mul(comp_log(y-ctr),box.b_l[0])
    for k = 1, p do
      particle.field = particle.field + comp_mul(comp_pow(y-ctr,k),
                                                  box.b_l[k])
--      c.printf("Box b_l[%i] = %f, %f", k, box.b_l[k]._1, box.b_l[k]._2)
    end

    -- Direct computation of influence from current box particles
    var num_near = box.num_part
    for j = 0, num_near do
      var id = box.part[j]
      if id ~= particle.id then
      var r : Vector2d = r_particles[box.part[j]].pos - particle.pos
      var kern : double = r:kfun()
      var field_new : Vector2d = r_particles[box.part[j]].q
                                 * cmath.pow(kern,2) * r
      particle.field = particle.field + field_new
      end
    end

    -- Direct computation of influence from neighbor box particles
    for j = 0, box.num_neighb do
      for k = 0, r_boxes[box.neighb[j]].num_part do
        var r : Vector2d = r_particles[r_boxes[box.neighb[j]].part[k]].pos
                           - particle.pos
        var kern : double = r:kfun()
        var field_new : Vector2d = r_particles[r_boxes[box.neighb[j]].part[k]].q
                                   * cmath.pow(kern,2) * r
        particle.field = particle.field + field_new

      end
    end

    particle.force = particle.q * particle.field
  end

end

task dump_forces(r_particles : region(Particle),
                 filename : int8[512],
                 token : uint32)
where
  reads(r_particles.field, r_particles.force)
do
  var f = c.fopen(filename,"w")
  for particle in r_particles do
    c.fprintf(f, "%18.16f %18.16f %18.16f %18.16f\n",
          particle.field._1, particle.field._2,
          particle.force._1, particle.force._2)
  end
  c.fclose(f)

  return 1
end

task dump_coeffs(r_boxes : region(ispace(int2d, box_per_dim_i2d), Box),
                 filename : int8[512],
                 p : uint32,
                 token : uint32)
where
  reads(r_boxes.a_k, r_boxes.b_l)
do
  var f = c.fopen(filename,"w")
  for box in r_boxes do
    for k = 0, p do
      c.fprintf(f, "%e %e ",
          box.a_k[k]._1, box.a_k[k]._2)
     end
     c.fprintf(f, "\n")
  end
  c.fclose(f)

  return 1
end

-- Main task
task toplevel()
  var ts_start_init = c.legion_get_current_time_in_micros()

  var config : FMMConfig
  config:initialize_from_command()

  -- Get force constant
  var force_const : double = get_force_const(config.model_type)
  c.printf("Force coefficient = %e\n", force_const)

  -- Create a region of particles
  var r_particles = region(ispace(ptr, config.num_particles), Particle)

  -- Allocate all the particles
  new(ptr(Particle, r_particles), config.num_particles)

  -- Initialize the particles from an input file
  initialize_particles(r_particles, config.input)

  var ts_stop_init = c.legion_get_current_time_in_micros()
  c.printf("Particle initialization took %.4f sec\n", (ts_stop_init - ts_start_init) * 1e-6)

  -- Setup variables for tree structure
  var ts_start_ilist = c.legion_get_current_time_in_micros()

  -- Calculate center and size of largest box
  var xmin : Vector2d = {1e32, 1e32}
  var xmax : Vector2d = {-1e32, -1e32}
  var qmax : double = 0.0
  for particle in r_particles do
    xmin._1 = min(xmin._1, particle.pos._1)
    xmin._2 = min(xmin._2, particle.pos._2)

    xmax._1 = max(xmax._1, particle.pos._1)
    xmax._2 = max(xmax._2, particle.pos._2)

    qmax = max(qmax, particle.q)
  end
  var ctr : Vector2d = 0.5*(xmin + xmax)
  var S : double = max(xmax._1 - xmin._1, xmax._2 - xmin._2)

  -- Calculate number of levels of refinement n = log_4(num_particles/3) to have
  -- ~3 particles per box in 2d, and calculate order p = log_2.12(epsilon)
  var num_lvl : uint32 = cmath.floor(cmath.log(config.num_particles/3)/cmath.log(4))
  var p : uint32 = cmath.ceil(cmath.log(config.epsilon/qmax)/cmath.log(2.12))

  p = sqrt(p*p)
  if p < 4 then p = 4
  elseif p > 20 then p = 20
  end

  c.printf("Number of mesh levels : %i\n", num_lvl)
  c.printf("Order p of multipole expansion : %i\n", p)

  -- Create a region of boxes
  var num_boxes : uint32 = 0
  var box_per_dim : uint32 = 0
  var box_per_dim_i2d : int2d = {0,0}
  for ibox = 0, num_lvl+1 do
    num_boxes = num_boxes + cmath.pow(4,ibox)
    box_per_dim = box_per_dim + cmath.pow(2,ibox)
  end
  box_per_dim_i2d = {box_per_dim,box_per_dim}
  c.printf("Num_boxes = %i\n", num_boxes)
  var r_boxes = region(ispace(int2d, box_per_dim_i2d), Box)

  -- Initialize interaction list tree structure
  init_Ilist_tree(r_particles, r_boxes, ctr, S, num_lvl, p)

  var ts_stop_ilist = c.legion_get_current_time_in_micros()
  c.printf("Interaction list setup took %.4f sec\n", (ts_stop_ilist - ts_start_ilist) * 1e-6)

  -- Upward pass
  var ts_start_up = c.legion_get_current_time_in_micros()
  var token : uint32
  token = M2M(r_particles, r_boxes, num_lvl, box_per_dim, p)
  var ts_stop_up = c.legion_get_current_time_in_micros()
  c.printf("Upward pass took %.4f sec\n", (ts_stop_up - ts_start_up) * 1e-6)

  -- Downward pass
  var ts_start_down = c.legion_get_current_time_in_micros()
  token = M2L(r_boxes, num_lvl, box_per_dim, p)
  c.printf("Evaluating forces\n")

  eval_forces(r_particles, r_boxes, num_lvl, box_per_dim, p)
  var ts_stop_down = c.legion_get_current_time_in_micros()
  c.printf("Downward pass took %.4f sec\n", (ts_stop_down - ts_start_down) * 1e-6)

--  token = dump_forces(r_particles, config.output, token)
  token = dump_coeffs(r_boxes, config.output, p, token)
end

regentlib.start(toplevel)
