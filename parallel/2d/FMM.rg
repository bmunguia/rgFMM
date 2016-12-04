
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
  phi      : double;   -- Potential
  boxes    : int2d[8]; -- Boxes containing particle
  id       : int32;    -- Particle ID
}

-- Fieldspace 'Box'
fspace Box {
  Parent     : int2d;      -- Parent box
  Children   : int2d[4];   -- Children boxes, up to 8 per box in 2d
  Ilist      : int2d[27]; -- Interaction list, up to 189 (6^2 - 3^2)
  neighb     : int2d[8];  -- Neighbors
  part       : int32[64];  -- Particles in box
  child_imin : int32;      -- Min i index of children
  child_imax : int32;      -- Max i index of children
  num_child  : int32;      -- Number of children
  num_ilist  : int32;      -- Number of interactions
  num_neighb : int32;      -- Number of neighbors
  num_part   : int32;      -- Number of particles in box
  ctr        : Vector2d;   -- Box center
  S          : double;     -- Box side length
  a_0        : double;     -- Coefficient in multipole expansion
  a_k        : Vector2d[24]; -- Coefficient in multipole expansion
  b_l        : Vector2d[24]; -- Coefficient in local expansion
}

terra skip_header(f : &c.FILE)
  var x0 : int8[512], y0 : int32, x1 : int8[512], y1 : int8[512]
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

terra Vector2d.metamethods.__div(v : Vector2d, c:double)
  return Vector2d { v._1 / c, v._2 / c }
end

terra well_separated(S : double, r : Vector2d)
  r._1 = cmath.round(sqrt(r._1 * r._1/(S*S)))
  r._2 = cmath.round(sqrt(r._2 * r._2/(S*S)))
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

terra comp_pow(vec : Vector2d, k : double)
  var theta : double = cmath.atan2(vec._2, vec._1)
  var z : double = sqrt(cmath.pow(vec._1,2)+cmath.pow(vec._1,2))
  z = cmath.pow(z,k)
  return Vector2d {z * cmath.cos(k*theta), z * cmath.sin(k*theta)}
end

terra comp_inv(vec : Vector2d)
  var tmp : Vector2d = Vector2d {vec._1, -vec._2}
  return tmp/(vec._1*tmp._1 - vec._2*tmp._2)
end

terra comp_log(vec : Vector2d)
  return Vector2d {cmath.log(sqrt(cmath.pow(vec._1,2) + cmath.pow(vec._2,2))),
                   cmath.atan2(vec._2,vec._1)}
end

terra factorial(n : int32)
  if n == 0 then return 1
  else
    var fact : int32 = 1
    for k = 1, n+1 do
      fact = fact * k
    end
    return fact
  end
end

terra choose(n : int32, k : int32)
  if k == 0 then return 1
  else
    var kf : int64 = 1
    var nf : int64 = 1
    for i = 2, k+1 do
      kf = kf * i
    end
    for i = n-k+1, n+1 do
      nf = nf * i
    end
    return nf/kf
  end
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

task create_tree(num_lvl     : int32,
                 r_particles : region(Particle),
                 r_boxes     : region(ispace(int2d, box_per_dim_i2d), Box),
                 ctr_0       : Vector2d,
                 S_0         : double)
where
  reads writes(r_boxes, r_particles)
do
  for particle in r_particles do
  var ctr = ctr_0
  var S = S_0
  
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
    var icell : int32 = int(x > ctr._1)
    var jcell : int32 = int(y > ctr._2)

    S = 0.5*S

    ctr._2 = ctr._2 + 0.5*S*(2.0*jcell-1.0)
    ctr._1 = ctr._1 + 0.5*S*(2.0*icell-1.0)

    -- Index of child node
    ibox = icell + 2*icell_old + offset
    jbox = jcell + 2*jcell_old + offset

    -- Check whether child node exists
    var child_exist : int8 = 0
    for j = 0, 4 do
      if r_boxes[icurr].Children[j] == [int2d]{ibox,jbox} then
        child_exist = 1
        break
      end
    end

    if child_exist == 0 then
      for j = 0, 4 do
        if r_boxes[icurr].Children[j] == [int2d]{-1,-1} then
        r_boxes[icurr].Children[j] = {ibox,jbox}
        r_boxes[icurr].num_child = r_boxes[icurr].num_child + 1
        r_boxes[r_boxes[icurr].Children[j]].Parent = icurr
        r_boxes[r_boxes[icurr].Children[j]].ctr = ctr
        r_boxes[r_boxes[icurr].Children[j]].S = S
	if ibox < r_boxes[icurr].child_imin then
	  r_boxes[icurr].child_imin = ibox
        end
	if jbox < r_boxes[icurr].child_imin then
	  r_boxes[icurr].child_imin = jbox
	end
	if ibox > r_boxes[icurr].child_imax then
	  r_boxes[icurr].child_imax = ibox
	end
	if jbox > r_boxes[icurr].child_imax then
	  r_boxes[icurr].child_imax = jbox
	end
	break
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
    icell_old = ibox - offset
    jcell_old = jbox - offset

    var box = r_boxes[{ibox,jbox}]

  end

  end

end

task create_neighb(r_boxes : region(ispace(int2d, box_per_dim_i2d), Box),
                   i       : int32,
		   j       : int32,
                   lvl     : int32)
where
  reads writes(r_boxes.num_neighb, r_boxes.neighb)
do
  var box_per_dim_lvl : int32 = cmath.pow(2,lvl)
  var offset : int32 = -1 + cmath.pow(2,lvl)

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

task create_Ilist(r_boxes : region(ispace(int2d), Box),
                  box_ind : int2d)
where
  reads(r_boxes),
  writes(r_boxes.num_ilist, r_boxes.Ilist)
do

  var S : double = r_boxes[box_ind].S
  var ctr : Vector2d = r_boxes[box_ind].ctr

  if r_boxes[box_ind].Parent ~= [int2d]{-1,-1} then
    var num_pn : int32 = r_boxes[r_boxes[box_ind].Parent].num_neighb
    var pn : int2d[8] = r_boxes[r_boxes[box_ind].Parent].neighb
    -- Loop over parent's neighbors
    for j = 0, num_pn do
      -- Loop over children
      for i = 0, r_boxes[pn[j]].num_child do
        var child : int2d = r_boxes[pn[j]].Children[i]
        var r : Vector2d = r_boxes[child].ctr - ctr
        if well_separated(S,r) then
          r_boxes[box_ind].num_ilist = r_boxes[box_ind].num_ilist + 1
          r_boxes[box_ind].Ilist[r_boxes[box_ind].num_ilist-1] = child
        end
      end
    end

  end


  for i = 0, r_boxes[box_ind].num_child do
    create_Ilist(r_boxes, r_boxes[box_ind].Children[i])
  end


end

task init_boxes(r_boxes : region(ispace(int2d, box_per_dim_i2d), Box),
                ctr     : Vector2d,
                S       : double,
                p       : int32)
where
  reads writes(r_boxes)
do

  -- Initialize boxes
  r_boxes[{0,0}].ctr = ctr
  r_boxes[{0,0}].S = S
  var child_init : int2d[4]
  var a_k_init : Vector2d[24]
  for i = 0, 4 do
    child_init[i] = int2d {-1, -1}
    a_k_init[i] = Vector2d {0.0, 0.0}
  end
  for i = 4, p do
    a_k_init[i] = Vector2d {0.0, 0.0}
  end
  fill(r_boxes.Children, child_init)
  fill(r_boxes.num_part, 0)
  fill(r_boxes.num_child, 0)
  fill(r_boxes.num_neighb, 0)
  fill(r_boxes.num_ilist, 0)
  fill(r_boxes.Parent, {-1, -1})
  fill(r_boxes.a_k, a_k_init)
  fill(r_boxes.b_l, a_k_init)
  fill(r_boxes.child_imin,1000)
  fill(r_boxes.child_imax,0)

end

task P2M(r_particles : region(Particle),
         r_boxes     : region(ispace(int2d, box_per_dim_i2d), Box),
         lvl         : int32,
         p           : int32)
where
  reads(r_particles.q, r_particles.pos, r_particles.boxes, r_boxes.ctr, r_boxes.a_k),
  writes(r_boxes.a_k)
do
  for particle in r_particles do
  
  r_boxes[particle.boxes[lvl]].a_k[0] = r_boxes[particle.boxes[lvl]].a_k[0]
                                       + Vector2d {particle.q, 0.0}

  var ds : Vector2d = particle.pos - r_boxes[particle.boxes[lvl]].ctr
  if ds._1 ~= 0.0 or ds._2 ~= 0.0 then
    var a_k : Vector2d = particle.q * Vector2d {1.0, 0.0}
    for k = 1, p do
      a_k = comp_mul(a_k,ds)
      r_boxes[particle.boxes[lvl]].a_k[k] = r_boxes[particle.boxes[lvl]].a_k[k]
                                            - (1/k) * a_k
    end
  end

  end

end

task outgoing_translation(r_boxes : region(ispace(int2d, box_per_dim_i2d), Box),
                          box_id  : int2d,
                          p       : int32)
where
  reads(r_boxes), writes(r_boxes.a_k)
do
  for k = 0, p do
    if k == 0 then
      for j = 0, r_boxes[box_id].num_child do
        r_boxes[box_id].a_k[k] = r_boxes[box_id].a_k[k]
                                 + r_boxes[r_boxes[box_id].Children[j]].a_k[k]
      end
    else
      for j = 0, r_boxes[box_id].num_child do
        var d : Vector2d = r_boxes[box_id].ctr
                           - r_boxes[r_boxes[box_id].Children[j]].ctr
        var a_sig : Vector2d = r_boxes[r_boxes[box_id].Children[j]].a_k[0]
	a_sig = comp_mul(comp_pow(d,k),a_sig)
        r_boxes[box_id].a_k[k] = r_boxes[box_id].a_k[k]
                                 - 1/k * a_sig
        for i = 1, k do
          var ch_pq : int32 = choose(k-1,i-1)
          a_sig = r_boxes[r_boxes[box_id].Children[j]].a_k[i]
	  a_sig = comp_mul(comp_pow(d,k-i),a_sig)
          r_boxes[box_id].a_k[k] = r_boxes[box_id].a_k[k]
                                   + ch_pq * a_sig
        end
      end
    end
  end
end

task incoming_translation(r_boxes : region(ispace(int2d, box_per_dim_i2d), Box),
                          box_id  : int2d,
                          p       : int32)
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
        var a_log_z : Vector2d = comp_log(Vector2d { -d._1,  -d._2})
	a_log_z = comp_mul(r_boxes[r_boxes[box_id].Ilist[j]].a_k[0],a_log_z)
        r_boxes[box_id].b_l[l] = r_boxes[box_id].b_l[l] + a_log_z
	var d_inv : Vector2d = comp_inv(d)
        for k = 1, p do
	  var b_mul : Vector2d = r_boxes[r_boxes[box_id].Ilist[j]].a_k[k]
          b_mul = comp_mul(comp_pow(d_inv,k),b_mul)
          r_boxes[box_id].b_l[l] = r_boxes[box_id].b_l[l]
                                   + cmath.pow(-1.0,k) * b_mul
        end
      end
      
    else
    
      for j = 0, r_boxes[box_id].num_ilist do
        var d : Vector2d =  r_boxes[r_boxes[box_id].Ilist[j]].ctr
                            - r_boxes[box_id].ctr

        var b_new : Vector2d = -1/l * r_boxes[r_boxes[box_id].Ilist[j]].a_k[0]
	var d_inv : Vector2d = comp_inv(d)
        for k = 1, p do
          var ch_pq : int32 = choose(l+k-1,k-1)
          b_new = b_new + ch_pq * cmath.pow(-1.0,k)
                  * comp_mul(r_boxes[r_boxes[box_id].Ilist[j]].a_k[k], comp_pow(d_inv,k))
        end
        r_boxes[box_id].b_l[l] = r_boxes[box_id].b_l[l]
                                 + comp_mul(b_new, comp_pow(d_inv,l))

      end
    end
--    c.printf("b_l = %e, %e\n", r_boxes[box_id].b_l[l]._1, r_boxes[box_id].b_l[l]._2)
  end

  end
  return 1
end

task L2L(r_boxes : region(ispace(int2d, box_per_dim_i2d), Box),
         box_id  : int2d,
         p       : int32)
where
  reads(r_boxes), writes(r_boxes.b_l)
do

  for j = 0, r_boxes[box_id].num_child do
    var box_child = r_boxes[box_id].Children[j]
    var d : Vector2d =  r_boxes[box_id].ctr - r_boxes[box_child].ctr
    for l = 0, p do
      for k = l, p do
        var ch_pq : int32 = choose(k,l)
        r_boxes[box_child].b_l[l] = r_boxes[box_child].b_l[l]
	                            + ch_pq * comp_mul(r_boxes[box_id].b_l[k],
				    comp_pow(Vector2d { -d._1,  -d._2}, k-l))
      end 
    end
  end

end

task M2M(r_boxes : region(ispace(int2d, box_per_dim_i2d), Box),
         box_min : int32,
         box_max : int32,
         p       : int32)
where
  reads(r_boxes), writes(r_boxes)
do
  -- Form multipole expansions about centers of all boxes at coarser mesh
  -- levels
  for ibox = box_min, box_max+1 do
    for jbox = box_min, box_max+1 do
      outgoing_translation(r_boxes, {ibox, jbox}, p)
    end
  end

  return 1

end

task M2L(r_boxes     : region(ispace(int2d, box_per_dim_i2d), Box),
         num_lvl     : int32,
         box_per_dim : int32,
         p           : int32)
where
  reads(r_boxes), writes(r_boxes.b_l)
do
  var box_min : int32 = 0
  var box_max : int32 = 0
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
  var token : int32 = 0
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
                 num_lvl     : int32,
                 box_per_dim : int32,
                 p           : int32)
where
  reads(r_particles, r_boxes),
  writes(r_particles.field, r_particles.force, r_particles.phi)
do

  for particle in r_particles do
    var box = r_boxes[particle.boxes[num_lvl]]
    var y = particle.pos
    var ctr = box.ctr
    var phi_vec = box.b_l[0]
    particle.field = particle.field + box.b_l[0]
    for k = 1, p do
      particle.field = particle.field + comp_mul(comp_pow(y-ctr,k),box.b_l[k])
      phi_vec = phi_vec + comp_mul(comp_pow(y-ctr,k),box.b_l[k])
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
      phi_vec = phi_vec + r_particles[box.part[j]].q * comp_log(r)
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
	phi_vec = phi_vec + r_particles[r_boxes[box.neighb[j]].part[k]].q * comp_log(r)

      end
    end

    particle.force = particle.q * particle.field
    particle.phi = phi_vec._1
  end

end

task dump_forces(r_particles : region(Particle),
                 filename : int8[512],
                 token : int32)
where
  reads(r_particles.field, r_particles.force, r_particles.phi)
do
  var f = c.fopen(filename,"w")
  for particle in r_particles do
    c.fprintf(f, "%18.16e\n", particle.phi)
  end
  c.fclose(f)

  return 1
end

task dump_coeffs(r_boxes : region(ispace(int2d, box_per_dim_i2d), Box),
                 filename : int8[512],
                 p : int32,
                 token : int32)
where
  reads(r_boxes.a_k, r_boxes.b_l)
do
  var f = c.fopen(filename,"w")
  for box in r_boxes do
    for k = 0, p do
      c.fprintf(f, "%e %e ",
          box.b_l[k]._1, box.b_l[k]._2)
          --box.a_k[k]._1, box.a_k[k]._2)
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
  var num_lvl : int32 = cmath.floor(cmath.log(config.num_particles/3)/cmath.log(4))
  var p : int32 = cmath.ceil(cmath.log(config.epsilon/qmax)/cmath.log(2.12))

  p = sqrt(p*p)
  if p < 4 then p = 4
  elseif p > 16 then p = 16
  end
  --num_lvl = 3

  c.printf("Number of mesh levels : %i\n", num_lvl)
  c.printf("Order p of multipole expansion : %i\n", p)

  -- Create a region of boxes
  var num_boxes : int32 = 0
  var box_per_dim : int32 = 0
  var box_per_dim_i2d : int2d = {0,0}
  for ilvl = 0, num_lvl+1 do
    num_boxes = num_boxes + cmath.pow(4,ilvl)
    box_per_dim = box_per_dim + cmath.pow(2,ilvl)
  end
  box_per_dim_i2d = {box_per_dim,box_per_dim}
  c.printf("Num_boxes = %i\n", num_boxes)
  var r_boxes = region(ispace(int2d, box_per_dim_i2d), Box)

  -- Create partitions
  var p_colors = ispace(int1d, config.parallelism)
  var p_particles = partition(equal, r_particles, p_colors)
  --var p_boxes = partition(equal, r_boxes, p_colors)

  -- Partition boxes based on levels
  var coloring = c.legion_domain_point_coloring_create()
  c.legion_domain_point_coloring_color_domain(coloring, [int1d](0), rect2d {{0,0},{0,0}})
  var i_min = 0
  var i_max = 0
  for ilvl = 1, num_lvl + 1 do
    i_min = i_min + cmath.pow(2,ilvl-1)
    i_max = i_max + cmath.pow(2,ilvl)
    c.legion_domain_point_coloring_color_domain(coloring, [int1d](ilvl),
                                                rect2d {{i_min,i_max},{i_min,i_max}})
  end
  var p_boxes_lvl = partition(disjoint, r_boxes, coloring, ispace(int1d, num_lvl+1))
  c.legion_domain_point_coloring_destroy(coloring)

  -- Initialize boxes
  init_boxes(r_boxes, ctr, S, p)

  -- Create FMM tree
  for color in p_colors do
    create_tree(num_lvl, p_particles[color], r_boxes, ctr, S)
  end
  
  -- Create neighbors list
  for lvl = 1, num_lvl+1 do
    var box_per_dim_lvl : int32 = cmath.pow(2,lvl)
    for i = 0, box_per_dim_lvl do
    for j = 0, box_per_dim_lvl do
      create_neighb(p_boxes_lvl[lvl], i, j, lvl)
    end
    end
  end

  -- Create interactions list
  create_Ilist(r_boxes, int2d {0,0})

  var ts_stop_ilist = c.legion_get_current_time_in_micros()
  c.printf("Interaction list setup took %.4f sec\n", (ts_stop_ilist - ts_start_ilist) * 1e-6)

  --
  -- Upward pass
  --
  var ts_start_up = c.legion_get_current_time_in_micros()
  
  -- Particle to multipole
  for color in p_colors do
    P2M(p_particles[color], p_boxes_lvl[num_lvl], num_lvl, p)
  end
  
  --Multipole to multipole
  var box_min : int32 = box_per_dim - cmath.pow(2,num_lvl)
  var box_max : int32 = box_per_dim - 1
  for lvl = num_lvl-1, -1, -1 do
    box_min = box_min - cmath.pow(2,lvl)
    box_max = box_max - cmath.pow(2,lvl+1)
    M2M(r_boxes, box_min, box_max, p)
  end

  var ts_stop_up = c.legion_get_current_time_in_micros()
  c.printf("Upward pass took %.4f sec\n", (ts_stop_up - ts_start_up) * 1e-6)

  --
  -- Downward pass
  --
  var ts_start_down = c.legion_get_current_time_in_micros()
  
  -- Multipole to local
  M2L(r_boxes, num_lvl, box_per_dim, p)

  -- Force evaluation
  c.printf("Evaluating forces\n")
  eval_forces(r_particles, r_boxes, num_lvl, box_per_dim, p)
  
  var ts_stop_down = c.legion_get_current_time_in_micros()
  c.printf("Downward pass took %.4f sec\n", (ts_stop_down - ts_start_down) * 1e-6)

  -- Dump output
  var token : int32
  for color in p_colors do
    token = dump_forces(p_particles[color], config.output, token)
  end
--  token = dump_coeffs(r_boxes, config.output, p, token)
end

regentlib.start(toplevel)
