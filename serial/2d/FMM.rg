
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
  id       : int64;    -- Particle ID
}

-- Fieldspace 'Box'
fspace Box {
  Parent     : int2d;      -- Parent box
  Children   : int2d[4];   -- Children boxes, up to 8 per box in 2d
  Ilist      : int2d[27]; -- Interaction list, up to 189 (6^2 - 3^2)
  neighb     : int2d[8];  -- Neighbors
  part       : int64[64];  -- Particles in box
  num_child  : int32;      -- Number of children
  num_ilist  : int32;      -- Number of interactions
  num_neighb : int32;      -- Number of neighbors
  num_part   : int32;      -- Number of particles in box
  ctr        : Vector2d;   -- Box center
  S          : double;     -- Box side length
  phi        : Vector2d;   -- p-term multipole expansion due to particles in box
  psi        : Vector2d;   -- p-term expansion due to interaction list
  psi_tilde  : Vector2d;   -- p-term local expansion due to particles in parent's
                           -- interaction list
  a_0        : double;     -- Coefficient in multipole expansion
  a_k        : Vector2d[4]; -- Coefficient in multipole expansion
}

terra skip_header(f : &c.FILE)
  var x0 : int8[512], y0 : uint64, x1 : int8[512], y1 : int8[512]
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

terra vec_pow(vec : Vector2d, k : uint64)
  var tmp1 : Vector2d = vec
  var tmp2 : Vector2d
  for i = 0, k do
    tmp2 = Vector2d {vec._1 * tmp1._1 - vec._2 * tmp1._2,
                     vec._1 * tmp1._2 + vec._2 * tmp1._1}
    tmp1 = tmp2
  end
  return tmp2
end

task initialize_particles(r_particles   : region(Particle),
                          filename      : int8[512])
where
  reads writes(r_particles)
do
  -- Set initial force and field to 0, and particle IDs
  var ipart : int64 = 0
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

task create_tree(num_lvl     : uint64,
                 r_particles : region(Particle),
                 r_boxes     : region(ispace(int2d, box_per_dim_i2d), Box),
                 particle    : ptr(Particle,r_particles),
                 ctr         : Vector2d,
                 S           : double)
where
  reads writes(r_particles, r_boxes)
do
  var x : double = particle.pos._1
  var y : double = particle.pos._2

  var box = r_boxes[{0,0}]
  particle.boxes[0] = {0,0}

  var icurr : int2d = {0, 0}
  var i : int32 = 0
  for lvl = 0, num_lvl do
    -- Index of child cell containing point
    var icell : uint32 = int(x > ctr._1)
    var jcell : uint32 = int(y > ctr._2)
    var ncell : uint32 = 2*icell + jcell

    S = 0.5*S

    ctr._2 = ctr._2 + 0.5*S*(2.0*jcell-1.0)
    ctr._1 = ctr._1 + 0.5*S*(2.0*icell-1.0)

    -- Index of child node
    var ibox : int32 = icell
    var jbox : int32 = jcell
    var nbox : int32 = ncell
    for j = 0, i do
      ibox = ibox + cmath.pow(2,j)
      jbox = jbox + cmath.pow(2,j)
    end

    -- Check whether child node exists
    if box.Children[ncell] == [int2d]{-1,-1} then
      box.Children[ncell] = {ibox,jbox}
      box.num_child = box.num_child + 1
      r_boxes[box.Children[ncell]].Parent = icurr
      r_boxes[box.Children[ncell]].ctr = ctr
      r_boxes[box.Children[ncell]].S = S
    end

    box.num_part = box.num_part + 1
    particle.boxes[lvl-1] = {ibox,jbox}
    box.part[box.num_part-1] = particle.id

    -- Update current box
    icurr = {ibox,jbox}
    i = i+1

    var box = r_boxes[{ibox,jbox}]

  end

end

task create_neighb(r_boxes : region(ispace(int2d, box_per_dim_i2d), Box),
                   num_lvl : uint64)
where
  reads writes(r_boxes)
do
  for lvl = 1,num_lvl do
  var box_per_dim_lvl : uint64 = cmath.pow(2,lvl)
  var offset : int64 = -1 + cmath.pow(2,lvl)
  for i = 0, box_per_dim_lvl do
  for j = 0, box_per_dim_lvl do
  -- Inner loop to find boxes within domain
    for ii = -1, 1, 2 do
    for jj = -1, 1, 2 do
      if (i+ii > -1) and (i+ii < box_per_dim_lvl) and (j+jj > -1) and
         (j+jj < box_per_dim_lvl) and (int2d {ii,jj} ~= int2d {0,0}) then
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
                  node    : Box)
where
  reads writes(r_boxes)
do
  var S : double = node.S
  var ctr : Vector2d = node.ctr

  if node.Parent ~= [int2d]{-1,-1} then
    var num_pn : int32 = r_boxes[node.Parent].num_neighb
    var pn : int2d[8] = r_boxes[node.Parent].neighb
    -- Loop over parent's neighbors
    for j = 0, num_pn do
      -- Loop over children
      for i = 0, 4 do
        var child : int2d = r_boxes[pn[j]].Children[i]
        var r : Vector2d = r_boxes[child].ctr - ctr
        if well_separated(S,r) then
          node.num_ilist = node.num_ilist + 1
          node.Ilist[node.num_ilist-1] = child
        end
      end
    end

  end

  for i = 0, 4 do
    if node.Children[i] ~= [int2d]{-1,-1} then
      create_Ilist(r_boxes, r_boxes[node.Children[i]])
    end
  end

end

task init_Ilist_tree(r_particles : region(Particle),
                     r_boxes     : region(ispace(int2d, box_per_dim_i2d), Box),
                     ctr         : Vector2d,
                     S           : double,
                     num_lvl     : uint64,
                     num_boxes   : uint64,
                     box_per_dim_i2d : int2d)
where
  reads writes(r_particles, r_boxes)
do
  var i : uint64 = 0

  -- Initialize boxes
  r_boxes[{0,0}].ctr = ctr
  r_boxes[{0,0}].S = S
  var child_init : int2d[4]
  var a_k_init : Vector2d[4]
  for i = 0, 4 do
    child_init[i] = {-1, -1}
    a_k_init[i] = {0.0, 0.0}
  end
  fill(r_boxes.Children, child_init)
  fill(r_boxes.num_part, 0)
  fill(r_boxes.num_child, 0)
  fill(r_boxes.num_neighb, 0)
  fill(r_boxes.Parent, {-1, -1})
  fill(r_boxes.phi, {0.0, 0.0})
  fill(r_boxes.psi, {0.0, 0.0})
  fill(r_boxes.psi_tilde, {0.0, 0.0})
  fill(r_boxes.a_k, a_k_init)
  -- Create FMM tree
  for particle in r_particles do
    create_tree(num_lvl, r_particles, r_boxes, particle, ctr, S)
  end

  -- Create neighbors list
  create_neighb(r_boxes, num_lvl)

  -- Create interactions list
  create_Ilist(r_boxes, r_boxes[{0,0}])

end

task P2M(r_particles   : region(Particle),
         r_boxes       : region(ispace(int2d, box_per_dim_i2d), Box),
         particle      : ptr(Particle,r_particles),
         lvl           : uint64,
         box_per_dim_lvl : uint64,
         p             : uint64)
where
  reads writes(r_particles, r_boxes)
do
  var ds : Vector2d = particle.pos - r_boxes[particle.boxes[lvl]].ctr

  if ds._1 ~= 0.0 or ds._2 ~= 0.0 then
    var a_k : Vector2d
    --phi = particle.q * phi
    for k = 0, p do
      a_k = -1/k * particle.q * vec_pow(ds,k)
      r_boxes[particle.boxes[lvl]].a_k[k] = r_boxes[particle.boxes[lvl]].a_k[k] + a_k
    end
    r_boxes[particle.boxes[lvl]].a_0 = r_boxes[particle.boxes[lvl]].a_0
                                       + particle.q
    return
  end

end

task M2M(r_particles : region(Particle),
         r_boxes     : region(ispace(int2d, box_per_dim_i2d), Box),
         num_lvl     : uint64,
         p           : uint64)
where
  reads writes(r_particles, r_boxes)
do
  -- Form multipole expansions of potential field due to particles in each box
  -- about the box center at the finest mesh level
  var lvl : uint64 = num_lvl
  var box_per_dim_lvl : uint64 = cmath.pow(2,lvl)
  for particle in r_particles do
    P2M(r_particles, r_boxes, particle, lvl, box_per_dim_lvl, p)
    --r_boxes[particle.boxes[lvl]].phi = r_boxes[particle.boxes[lvl]].phi + phi
  end

  -- Form multipole expansions about centers of all boxes at coarser mesh
  -- levels
  for lvl = num_lvl-1, -1, 0 do
    for box in r_boxes[particle.boxes[lvl]] do

    end
  end

end

task dump_forces(r_particles : region(Particle),
                 filename : int8[512])
where
  reads(r_particles)
do
  var f = c.fopen(filename,"w")
  for particle in r_particles do
    c.fprintf(f, "%f %f %f\n",
          particle.field, particle.force._1, particle.force._2)
  end
  c.fclose(f)
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
  for particle in r_particles do
    xmin._1 = min(xmin._1, particle.pos._1)
    xmin._2 = min(xmin._2, particle.pos._2)

    xmax._1 = max(xmax._1, particle.pos._1)
    xmax._2 = max(xmax._2, particle.pos._2)
  end
  var ctr : Vector2d = 0.5*(xmin + xmax)
  var S : double = max(xmax._1 - xmin._1, xmax._2 - xmin._2)

  -- Calculate number of levels of refinement n = log_4(num_particles/3) to have
  -- ~3 particles per box in 2d, and calculate order p = log_2(epsilon)
  var num_lvl : uint64 = cmath.floor(cmath.log(config.num_particles/3)/cmath.log(4))
  var p : uint64 = cmath.ceil(cmath.log(1/config.epsilon)/cmath.log(2))
  if p < 1 then p = 1 end

  c.printf("Number of mesh levels : %i\n", num_lvl)
  c.printf("Order p of multipole expansion : %i\n", p)

  -- Create a region of boxes
  var num_boxes : uint64 = 0
  var box_per : uint64 = 0
  var box_per_dim_i2d : int2d = {0,0}
  for ibox = 0, num_lvl+1 do
    num_boxes = num_boxes + cmath.pow(4,ibox)
    box_per = box_per + cmath.pow(2,ibox)
  end
  box_per_dim_i2d = {box_per,box_per}
  c.printf("Num_boxes = %i\n", num_boxes)
  var r_boxes = region(ispace(int2d, box_per_dim_i2d), Box)

  -- Initialize interaction list tree structure
  init_Ilist_tree(r_particles, r_boxes, ctr, S, num_lvl, num_boxes, box_per_dim_i2d)

  var ts_stop_ilist = c.legion_get_current_time_in_micros()
  c.printf("Interaction list setup took %.4f sec\n", (ts_stop_ilist - ts_start_ilist) * 1e-6)

  -- Upward pass
  var ts_start_up = c.legion_get_current_time_in_micros()
  M2M(r_particles, r_boxes, num_lvl, p)
  var ts_stop_up = c.legion_get_current_time_in_micros()
  c.printf("Upward pass took %.4f sec\n", (ts_stop_up - ts_start_up) * 1e-6)


  -- Downward pass

--  c.printf("Force calculation took %.4f sec\n", (ts_stop_forces - ts_start_forces) * 1e-6)

  dump_forces(r_particles, config.output)
end

regentlib.start(toplevel)
