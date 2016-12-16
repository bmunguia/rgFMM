
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

struct minmax {
  _1 : double;  --xmin
  _2 : double;  --xmax
  _3 : double;  --ymin
  _4 : double;  --ymax
  _5 : double;  --qmax
}

-- Fieldspace 'Box'
fspace Box {
  Parent     : int2d;      -- Parent box
  Children   : int2d[4];   -- Children boxes, up to 8 per box in 2d
  Ilist      : int2d[27];  -- Interaction list, up to 189 (6^2 - 3^2)
  neighb     : int2d[8];   -- Neighbors
  part       : ptr[24];  -- Particles in box
  child_imin : int64;      -- Min i index of children
  child_imax : int64;      -- Max i index of children
  child_jmin : int64;      -- Min j index of children
  child_jmax : int64;      -- Max j index of children
  ilist_imin : int64;      -- Min i index of Ilist
  ilist_imax : int64;      -- Max i index of Ilist
  ilist_jmin : int64;      -- Min j index of Ilist
  ilist_jmax : int64;      -- Max j index of Ilist
  num_child  : int64;      -- Number of children
  num_ilist  : int64;      -- Number of interactions
  num_neighb : int64;      -- Number of neighbors
  num_part   : int64;      -- Number of particles in box
  ctr        : Vector2d;   -- Box center
  S          : double;     -- Box side length
  a_k        : Vector2d[20]; -- Coefficient in multipole expansion
  b_l        : Vector2d[20]; -- Coefficient in local expansion
}

-- Fieldspace 'Particle'
fspace Particle(r: region(ispace(int2d), Box)) {
  q        : double;   -- Particle source strength
  pos      : Vector2d; -- Position vector
  vel      : Vector2d; -- Velocity vector
  force    : Vector2d; -- Force vector
  field    : Vector2d; -- Field strength
  phi      : double;   -- Potential
  boxes    : int2d(Box,r); -- Leaf box containing particle
  id       : ptr;    -- Particle ID
}

task factorize(parallelism : int) : int2d
  var limit = [int](cmath.sqrt([double](parallelism)))
  var size_x = 1
  var size_y = parallelism
  for i = 1, limit + 1 do
    if parallelism % i == 0 then
      size_x, size_y = i, parallelism / i
      if size_x > size_y then
        size_x, size_y = size_y, size_x
      end
    end
  end
  return int2d { size_x, size_y }
end

terra wait_for(x : int) return 1 end

terra skip_header(f : &c.FILE)
  var x0 : int8[512], y0 : int64, x1 : int8[512], y1 : int8[512]
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
  return c.fscanf(f, "%s %lf %lf %lf %lf %lf\n", &x, &source[0],
                  &source[1], &source[2]) == 4
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
  var r1 : double = cmath.round(sqrt(r._1 * r._1/(S*S)))
  var r2 : double = cmath.round(sqrt(r._2 * r._2/(S*S)))
  if r1 > 1.0 or r2 > 1.0 then
    return true
  else
    return false
  end
end

terra comp_mul(vec1 : Vector2d, vec2 : Vector2d)
  return Vector2d {vec1._1 * vec2._1 - vec1._2 * vec2._2,
                   vec1._2 * vec2._1 + vec1._1 * vec2._2}
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

terra factorial(n : int64)
  if n == 0 then return 1
  else
    var fact : int64 = 1
    for k = 1, n+1 do
      fact = fact * k
    end
    return fact
  end
end

terra choose(n : int64, k : int64)
  if k == 0 then return 1
  elseif k == n then return 1
  else
    var nf : uint64 = factorial(n)
    var kf : uint64 = factorial(k)*factorial(n-k)
    return nf/kf
  end
end

task initialize_particles(r_boxes     : region(ispace(int2d), Box),
                          r_particles : region(Particle(r_boxes)),
                          filename    : int8[512])
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
    --particle.vel = Vector2d {source[3], source[4]}
  end

  c.fclose(f)

end

task get_minmax(r_boxes     : region(ispace(int2d), Box),
                r_particles : region(Particle(r_boxes)),
                minmax_old  : minmax)
where
  reads(r_particles.{pos, q})
do
  for particle in r_particles do
    minmax_old._1 = min(minmax_old._1, particle.pos._1)
    minmax_old._2 = max(minmax_old._2, particle.pos._1)
    minmax_old._3 = min(minmax_old._3, particle.pos._2)
    minmax_old._4 = max(minmax_old._4, particle.pos._2)
    minmax_old._5 = max(minmax_old._5, sqrt(particle.q*particle.q))
  end
  return minmax_old
      
end

task create_tree(num_lvl      : int64,
                 r_boxes      : region(ispace(int2d), Box),
		 r_boxes_leaf : region(ispace(int2d), Box),
		 r_particles  : region(Particle(r_boxes_leaf)),
                 ctr_0        : Vector2d,
                 S_0          : double)
where
  reads writes(r_boxes, r_particles)
do
  for particle in r_particles do
  var ctr = ctr_0
  var S = S_0

  var x : double = particle.pos._1
  var y : double = particle.pos._2

--  particle.boxes[0] = {0,0}

  var icurr : int2d = {0, 0}
  var ibox : int64 = 0
  var jbox : int64 = 0
  var icell_old : int64 = 0
  var jcell_old : int64 = 0
  var i : int64 = 1
  var offset : int64 = -1
  for lvl = 0, num_lvl do
    -- Index of child cell containing point
    var offset : int64 = -1 + cmath.pow(2,lvl+1)
    var icell : int64 = int(x > ctr._1)
    var jcell : int64 = int(y > ctr._2)

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
	if jbox < r_boxes[icurr].child_jmin then
	  r_boxes[icurr].child_jmin = jbox
	end
	if ibox > r_boxes[icurr].child_imax then
	  r_boxes[icurr].child_imax = ibox
	end
	if jbox > r_boxes[icurr].child_imax then
	  r_boxes[icurr].child_jmax = jbox
	end
	break
        end
      end
    end

    if lvl == num_lvl-1 then
      r_boxes[{ibox,jbox}].num_part = r_boxes[{ibox,jbox}].num_part + 1
      r_boxes[{ibox,jbox}].part[r_boxes[{ibox,jbox}].num_part-1] =
                           unsafe_cast(ptr(Particle(r_boxes_leaf),r_particles),particle.id)
      particle.boxes = unsafe_cast(int2d(Box,r_boxes_leaf), {ibox,jbox})
    end

    -- Update current box
    icurr = {ibox,jbox}
    i = i+1
    icell_old = ibox - offset
    jcell_old = jbox - offset

    var box = r_boxes[{ibox,jbox}]

  end

  end

end

task create_neighb(r_boxes : region(ispace(int2d), Box),
                   lvl     : int64)
where
  reads writes(r_boxes)
do
  var box_per_dim_lvl : int64 = cmath.pow(2,lvl)
  var offset : int64 = -1 + cmath.pow(2,lvl)
  var r_bounds = r_boxes.bounds

  for i = 0, box_per_dim_lvl do
  for j = 0, box_per_dim_lvl do

  -- Inner loop to find boxes within domain
  for ii = -1, 2 do
  for jj = -1, 2 do
    if (i+ii > -1) and (i+ii < box_per_dim_lvl) and (j+jj > -1) and
       (j+jj < box_per_dim_lvl) and (ii ~= 0 or jj ~= 0) then
         var num_neighb : int64 = r_boxes[{i+offset, j+offset}].num_neighb + 1
         r_boxes[{i+offset, j+offset}].num_neighb = num_neighb
         r_boxes[{i+offset, j+offset}].neighb[num_neighb-1] = {i+ii+offset, j+jj+offset}
    end
  end
  end

  end
  end

end

task create_Ilist(r_boxes : region(ispace(int2d), Box),
                  box_ind : int2d)
where
  reads(r_boxes),
  writes(r_boxes.num_ilist, r_boxes.Ilist, r_boxes.ilist_imin,
         r_boxes.ilist_imax, r_boxes.ilist_jmin, r_boxes.ilist_jmax)
do

  var S : double = r_boxes[box_ind].S
  var ctr : Vector2d = r_boxes[box_ind].ctr

  if r_boxes[box_ind].Parent ~= [int2d]{-1,-1} then
    var num_pn : int64 = r_boxes[r_boxes[box_ind].Parent].num_neighb
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
      -- Get min and max indices for interaction list partitioning
      if r_boxes[pn[j]].child_imin < r_boxes[box_ind].ilist_imin then
	r_boxes[box_ind].ilist_imin = r_boxes[pn[j]].child_imin
      end
      if r_boxes[pn[j]].child_imax > r_boxes[box_ind].ilist_imax then
	r_boxes[box_ind].ilist_imax = r_boxes[pn[j]].child_imax
      end
      if r_boxes[pn[j]].child_jmin < r_boxes[box_ind].ilist_jmin then
	r_boxes[box_ind].ilist_jmin = r_boxes[pn[j]].child_jmin
      end
      if r_boxes[pn[j]].child_jmax > r_boxes[box_ind].ilist_jmax then
	r_boxes[box_ind].ilist_jmax = r_boxes[pn[j]].child_jmax
      end

    end

  end


  for i = 0, r_boxes[box_ind].num_child do
    create_Ilist(r_boxes, r_boxes[box_ind].Children[i])
  end


end

task create_Ilist_new(r_boxes     : region(ispace(int2d), Box),
                      r_boxes_par : region(ispace(int2d), Box))
where
  reads(r_boxes, r_boxes_par),
  writes(r_boxes.num_ilist, r_boxes.Ilist, r_boxes.ilist_imin,
         r_boxes.ilist_imax, r_boxes.ilist_jmin, r_boxes.ilist_jmax)
do
  for box in r_boxes do
  
  var S : double = box.S
  var ctr : Vector2d = box.ctr

  if box.Parent ~= [int2d]{-1,-1} then
    var num_pn : int64 = r_boxes_par[box.Parent].num_neighb
    var pn : int2d[8] = r_boxes_par[box.Parent].neighb
    -- Loop over parent's neighbors
    for j = 0, num_pn do
      -- Loop over children
      for i = 0, r_boxes_par[pn[j]].num_child do
        var child : int2d = r_boxes_par[pn[j]].Children[i]
        var r : Vector2d = r_boxes[child].ctr - ctr
        if well_separated(S,r) then
          box.num_ilist = box.num_ilist + 1
          box.Ilist[box.num_ilist-1] = child
        end
      end
      -- Get min and max indices for interaction list partitioning
      if r_boxes_par[pn[j]].child_imin < box.ilist_imin then
	box.ilist_imin = r_boxes_par[pn[j]].child_imin
      end
      if r_boxes_par[pn[j]].child_imax > box.ilist_imax then
	box.ilist_imax = r_boxes_par[pn[j]].child_imax
      end
      if r_boxes_par[pn[j]].child_jmin < box.ilist_jmin then
	box.ilist_jmin = r_boxes_par[pn[j]].child_jmin
      end
      if r_boxes_par[pn[j]].child_jmax > box.ilist_jmax then
	box.ilist_jmax = r_boxes_par[pn[j]].child_jmax
      end

    end

  end
  
  end
end

task init_boxes(r_boxes : region(ispace(int2d), Box),
                ctr     : Vector2d,
                S       : double,
                p       : int64)
where
  reads writes(r_boxes)
do

  -- Initialize boxes
  r_boxes[{0,0}].ctr = ctr
  r_boxes[{0,0}].S = S
  var child_init : int2d[4]
  var a_k_init : Vector2d[20]
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
  fill(r_boxes.child_imin,1e6)
  fill(r_boxes.child_imax,-1)
  fill(r_boxes.child_jmin,1e6)
  fill(r_boxes.child_jmax,-1)
  fill(r_boxes.ilist_imin,1e6)
  fill(r_boxes.ilist_imax,-1)
  fill(r_boxes.ilist_jmin,1e6)
  fill(r_boxes.ilist_jmax,-1)

end

task P2M(r_boxes     : region(ispace(int2d), Box),
	 r_particles : region(Particle(r_boxes)),
         lvl         : int64,
         p           : int64)
where
  reads(r_particles.{q, pos, boxes}, r_boxes.{ctr, a_k}),
  writes(r_boxes.a_k)
do
  for particle in r_particles do

  r_boxes[particle.boxes].a_k[0] = r_boxes[particle.boxes].a_k[0]
                                       + Vector2d {particle.q, 0.0}

  var ds : Vector2d = particle.pos - r_boxes[particle.boxes].ctr
  if ds._1 ~= 0.0 or ds._2 ~= 0.0 then
    for k = 1, p do
      r_boxes[particle.boxes].a_k[k] += - (1/k) * particle.q * comp_pow(ds,k)
    end
  end

  end
  return 1

end

task M2M(r_boxes       : region(ispace(int2d), Box),
         r_boxes_child : region(ispace(int2d), Box),
         p             : int64)
where
  reads(r_boxes.{num_child, a_k, ctr, Children}, r_boxes_child.{a_k, ctr}), writes(r_boxes.a_k)
do
  
  for box in r_boxes do
  if box.num_child > 0 then

  for k = 0, p do
    if k == 0 then
      for j = 0, box.num_child do
        var a_sig = r_boxes_child[box.Children[j]].a_k[k]
        box.a_k[k] += a_sig
      end
    else
      for j = 0, box.num_child do
        var d : Vector2d = r_boxes_child[box.Children[j]].ctr - box.ctr
        var a_sig : Vector2d = r_boxes_child[box.Children[j]].a_k[0]
	a_sig = comp_mul(comp_pow(d,k),a_sig)
        box.a_k[k] += - 1/k * a_sig
        for i = 1, k+1 do
          var ch_pq : int64 = choose(k-1,i-1)
          a_sig = r_boxes_child[box.Children[j]].a_k[i]
	  a_sig = comp_mul(comp_pow(d,k-i),a_sig)
          box.a_k[k] += ch_pq * a_sig
        end
      end
    end
  end

  end
  end
  return 1

end

task M2L(r_boxes       : region(ispace(int2d), Box),
         r_boxes_ilist : region(ispace(int2d), Box),
         p             : int64)
where
  reads(r_boxes.{num_ilist, Ilist, ctr, b_l}, r_boxes_ilist.{ctr, a_k}), writes(r_boxes.b_l)
do

  for box in r_boxes do
  if box.num_ilist > 0 then

  for l = 0, p do
    if l == 0 then
      var b_l = Vector2d {0.0, 0.0}
      for j = 0, box.num_ilist do
        var d : Vector2d =  r_boxes_ilist[box.Ilist[j]].ctr - box.ctr
        var a_log_z : Vector2d = comp_log(Vector2d { -d._1,  -d._2})
	a_log_z = comp_mul(r_boxes_ilist[box.Ilist[j]].a_k[0],a_log_z)
        b_l = b_l + a_log_z
	var d_inv : Vector2d = comp_inv(d)
	var b_new = Vector2d {0.0, 0.0}
        for k = 1, p do
	  var b_mul : Vector2d = r_boxes_ilist[box.Ilist[j]].a_k[k]
          b_mul = comp_mul(comp_pow(d_inv,k),b_mul)
          b_l = b_l +  cmath.pow(-1.0,k) * b_mul
        end
      end
      box.b_l[l] = b_l

    else
      var b_l = Vector2d {0.0, 0.0}
      for j = 0, box.num_ilist do
        var b_new = Vector2d {0.0, 0.0}
        var d : Vector2d =  r_boxes_ilist[box.Ilist[j]].ctr - box.ctr
        b_new =  b_new - 1/l * r_boxes_ilist[box.Ilist[j]].a_k[0]
	var d_inv : Vector2d = comp_inv(d)
        for k = 1, p do
          var ch_pq : int64 = choose(l+k-1,k-1)
	  if ch_pq < 0 then
 	    c.printf("ch_pq = %i\n", ch_pq)
	  end
          b_new = b_new + ch_pq * cmath.pow(-1.0,k)
                  * comp_mul(r_boxes_ilist[box.Ilist[j]].a_k[k], comp_pow(d_inv,k))
        end
	b_l = b_l + comp_mul(b_new, comp_pow(d_inv,l))
      end
      box.b_l[l] = b_l
    end
  end

  end
  end
  return 1
end

task L2L(r_boxes       : region(ispace(int2d), Box),
         r_boxes_child : region(ispace(int2d), Box),
         p             : int64)
where
  reads(r_boxes.{ctr, b_l, num_child, Children}, r_boxes_child.{ctr, b_l}), writes(r_boxes_child.b_l)
do
  for box in r_boxes do
  for j = 0, box.num_child do
    var d : Vector2d =  box.ctr - r_boxes_child[box.Children[j]].ctr
    for l = 0, p do
      for k = l, p do
        var ch_pq : int64 = choose(k,l)
        r_boxes_child[box.Children[j]].b_l[l] += ch_pq * comp_mul(box.b_l[k],
				           comp_pow(Vector2d { -d._1,  -d._2}, k-l))
      end
    end
  end
  end

end

task eval_forces(r_boxes            : region(ispace(int2d), Box),
                 r_boxes_neighb     : region(ispace(int2d), Box),
		 r_particles        : region(Particle(r_boxes)),
		 r_particles_neighb : region(Particle(r_boxes_neighb)),
                 num_lvl            : int64,
                 box_per_dim        : int64,
                 p                  : int64)
where
  reads(r_boxes.{num_part, part, b_l, ctr, num_neighb, neighb}, r_boxes_neighb.{ctr, num_part, part}, r_particles_neighb.{pos, q, id}),
  writes(r_particles.phi)
do

  for box in r_boxes do
  for j = 0, box.num_part do
    var y = r_particles_neighb[box.part[j]].pos
    var ctr = box.ctr
    var phi_vec = box.b_l[0]
--    particle.field = particle.field + box.b_l[0]
    for k = 1, p do
--      particle.field = particle.field + comp_mul(comp_pow(y-ctr,k),box.b_l[k])
      phi_vec = phi_vec + comp_mul(comp_pow(y-ctr,k),box.b_l[k])
    end

    -- Direct computation of influence from current box particles
    var num_near = box.num_part
    for i = 0, num_near do
--      if r_particles_neighb[box.part[i]].id ~= r_particles_neighb[box.part[j]].id then
      var r : Vector2d = y - r_particles_neighb[box.part[i]].pos
--      var kern : double = r:kfun()
--      var field_new : Vector2d = r_particles_neighb[box.part[i]].q
--                                 * cmath.pow(kern,2) * r
      if r._1 ~= 0.0 or r._2 ~= 0.0 then
--      particle.field = particle.field + field_new
        phi_vec = phi_vec + r_particles_neighb[box.part[i]].q * comp_log(r)
      end
--      end
    end

    -- Direct computation of influence from neighbor box particles
    for i = 0, box.num_neighb do
      for k = 0, r_boxes_neighb[box.neighb[i]].num_part do
        var r : Vector2d = y - r_particles_neighb[r_boxes_neighb[box.neighb[i]].part[k]].pos
--        var kern : double = r:kfun()
--        var field_new : Vector2d = r_particles_neighb[r_boxes_neighb[box.neighb[i]].part[k]].q
--                                   * cmath.pow(kern,2) * r
--        particle.field = particle.field + field_new
--        if r._1 ~= 0.0 or r._2 ~= 0.0 then
	  phi_vec = phi_vec + r_particles_neighb[r_boxes_neighb[box.neighb[i]].part[k]].q * comp_log(r)
--	end
      end
    end

--    r_particles[box.part[j]].field = particle.field
--    r_particles[box.part[j]].force = particle.q * particle.field
    r_particles[box.part[j]].phi = phi_vec._1
  end

  end

  return 1

end

task dump_forces(r_boxes     : region(ispace(int2d), Box),
                 r_particles : region(Particle(r_boxes)),
                 filename : int8[512],
                 token : int32)
where
  reads(r_particles.phi, r_particles.id, r_particles.pos)
do
  var f = c.fopen(filename,"a")
  for particle in r_particles do
    c.fprintf(f, "%i, %18.16f, %18.16f, %18.16f\n", particle.id, particle.phi, particle.pos._1, particle.pos._2)
  end
  c.fclose(f)

  return 1
end

task dump_coeffs(r_boxes : region(ispace(int2d, box_per_dim_i2d), Box),
                 filename : int8[512],
                 p : int64,
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

  -- Calculate number of levels of refinement n = log_4(num_particles/3) to have
  -- ~3 particles per box in 2d
  var num_lvl : int64 = cmath.floor(cmath.log(config.num_particles/3)/cmath.log(4))

  c.printf("Number of mesh levels : %i\n", num_lvl)

  -- Create a region of boxes
  var num_boxes : int64 = 0
  var box_per_dim : int64 = 0
  var box_per_dim_i2d : int2d = {0,0}
  for ilvl = 0, num_lvl+1 do
    num_boxes = num_boxes + cmath.pow(4,ilvl)
    box_per_dim = box_per_dim + cmath.pow(2,ilvl)
  end
  box_per_dim_i2d = {box_per_dim,box_per_dim}
  c.printf("Num_boxes = %i\n", num_boxes)

  var r_boxes = region(ispace(int2d, box_per_dim_i2d), Box)

  -- Partition boxes based on levels
  var coloring = c.legion_domain_point_coloring_create()
  c.legion_domain_point_coloring_color_domain(coloring, [int1d](0), rect2d {{0,0},{0,0}})
  var box_min : int64[12]
  var box_max : int64[12]
  var i_min = 0
  var i_max = 0
  box_min[0] = i_min
  box_max[0] = i_max
  for ilvl = 1, num_lvl + 1 do
    i_min = i_min + cmath.pow(2,ilvl-1)
    i_max = i_max + cmath.pow(2,ilvl)
    box_min[ilvl] = i_min
    box_max[ilvl] = i_max
    c.legion_domain_point_coloring_color_domain(coloring, [int1d](ilvl),
                                                rect2d {{i_min,i_min},{i_max,i_max}})
  end
  var p_boxes_lvl = partition(disjoint, r_boxes, coloring, ispace(int1d, num_lvl+1))
  c.legion_domain_point_coloring_destroy(coloring)

  var r_boxes_leaf = p_boxes_lvl[num_lvl]

  -- Create a region of particles
  var r_particles = region(ispace(ptr, config.num_particles), Particle(wild))

  -- Allocate all the particles
  new(ptr(Particle(r_boxes_leaf), r_particles), config.num_particles)

  -- Create partitions
  var p_colors = ispace(int1d, config.parallelism)
--  var p_boxes_part = partition(equal, r_boxes_leaf, p_colors)
--  var p_particles = preimage(r_particles, p_boxes_part, r_particles.boxes)
  var p_particles = partition(equal, r_particles, p_colors)

  -- Initialize the particles from an input file
  var xmin : Vector2d = {1e32, 1e32}
  var xmax : Vector2d = {-1e32, -1e32}
  var qmax : double = 0.0
  --for color in p_colors do
  --  initialize_particles(p_boxes_part[color], p_particles[color], config.input)
  --end
  initialize_particles(r_boxes_leaf, r_particles, config.input)
  var root_minmax : minmax = {xmin._1, xmax._1, xmin._2, xmax._2, qmax}
  var minmax_new : minmax = root_minmax
  for color in p_colors do
--    for particle in p_particles[color] do
      --xmin._1 = min(xmin._1, particle.pos._1)
      --xmin._2 = min(xmin._2, particle.pos._2)
      --xmax._1 = max(xmax._1, particle.pos._1)
      --xmax._2 = max(xmax._2, particle.pos._2)
      --qmax = max(qmax, sqrt(particle.q*particle.q))
      minmax_new = get_minmax(r_boxes_leaf, p_particles[color], root_minmax)
      root_minmax = minmax_new
--    end
  end
  xmin = Vector2d{root_minmax._1, root_minmax._3}
  xmax = Vector2d{root_minmax._2, root_minmax._4}
  qmax = root_minmax._5
  var ctr : Vector2d = 0.5*(xmin + xmax)
  var S : double = max(xmax._1 - xmin._1, xmax._2 - xmin._2)

  var ts_stop_init = c.legion_get_current_time_in_micros()
  c.printf("Particle initialization took %.4f sec\n", (ts_stop_init - ts_start_init) * 1e-6)

  -- Setup variables for tree structure
  var ts_start_ilist = c.legion_get_current_time_in_micros()
  
  -- Calculate order p = log_2.12(epsilon)
  var p : int64 = cmath.ceil(cmath.log(config.epsilon/qmax)/cmath.log(2.12))
  p = sqrt(p*p)
  if p < 4 then p = 4
  elseif p > 20 then p = 20
  end
  c.printf("Order p of multipole expansion : %i\n", p)

  -- Initialize boxes
  init_boxes(r_boxes, ctr, S, p)

  -- Create FMM tree
  for color in p_colors do
    create_tree(num_lvl, r_boxes, r_boxes_leaf, p_particles[color], ctr, S)
  end

--  var p_boxes_part = image(r_boxes_leaf, p_particles, r_particles.boxes)

  -- Create neighbors list
  for lvl = 1, num_lvl+1 do
      create_neighb(p_boxes_lvl[lvl], lvl)
  end

  -- Create interactions list
--  create_Ilist(r_boxes, int2d {0,0})
  for lvl = 2, num_lvl+1 do
    create_Ilist_new(p_boxes_lvl[lvl], p_boxes_lvl[lvl-1])
  end

  var ts_stop_ilist = c.legion_get_current_time_in_micros()
  c.printf("Interaction list setup took %.4f sec\n", (ts_stop_ilist - ts_start_ilist) * 1e-6)

  -- Add partitions which include box children, parents, neighbors, and Ilist
  var c_self = c.legion_domain_point_coloring_create()
  var c_ilist = c.legion_domain_point_coloring_create()

  var icolor : int64 = 0
  var ileaf  : int64 = 0
  for ilvl = 1, num_lvl + 1 do
    i_min = box_min[ilvl]
    i_max = box_max[ilvl]
    var r = p_boxes_lvl[ilvl]
    if ilvl  < 5 then
      var box_min = r[{i_min,i_min}]
      var box_max = r[{i_max,i_max}]
      c.legion_domain_point_coloring_color_domain(c_self, [int1d](icolor),
            rect2d{{i_min,i_min},{i_max,i_max}})
      c.legion_domain_point_coloring_color_domain(c_ilist, [int1d](icolor),
            rect2d{{i_min,i_min},{i_max,i_max}})
      icolor += 1
    else
      var j_init : int64 = i_min
      var j_inc : int64 = cmath.pow(2,5)
      for j_min = j_init, i_max, j_inc do
        for ii_min = j_init, i_max, j_inc do
          var j_max = [int64](j_min+j_inc)
	  var ii_max = [int64](ii_min+j_inc)
	  if j_max > i_max then
	    j_max = i_max
	  end
	  if ii_max > i_max then
	    ii_max = i_max
	  end
          c.legion_domain_point_coloring_color_domain(c_self, [int1d](icolor),
               rect2d{{ii_min,j_min},{ii_max,j_max}})
          var ilist_imin = ii_min-3
	  var ilist_jmin = j_min-3
	  var ilist_imax = ii_max+3
	  var ilist_jmax = j_max+3
	  if ilist_imin < i_min then
	    ilist_imin = i_min
	  end
	  if ilist_imax > i_max then
	    ilist_imax = i_max
	  end
	  if ilist_jmin < i_min then
	    ilist_jmin = i_min
	  end
	  if ilist_jmax > i_max then
	    ilist_jmax = i_max
	  end
	  c.legion_domain_point_coloring_color_domain(c_ilist, [int1d](icolor),
               rect2d{{ilist_imin,ilist_jmin},{ilist_imax,ilist_jmax}})
          icolor += 1
	end
      end
    end
  end
  var p_boxes_self = partition(disjoint, r_boxes, c_self, ispace(int1d, icolor))
  var p_boxes_ilist = partition(aliased, r_boxes, c_ilist, ispace(int1d, icolor))
  c.legion_domain_point_coloring_destroy(c_self)
  c.legion_domain_point_coloring_destroy(c_ilist)

  -- Partition particles based on box at leaf level, and by neighbors
  var p_colors_boxes = ispace(int2d, factorize(config.parallelism))
  var p_boxes_leaf = partition(equal, r_boxes_leaf, p_colors_boxes)
  var c_neighb = c.legion_domain_point_coloring_create()
  for color in p_colors_boxes do
    var bounds_lo = p_boxes_leaf[color].bounds.lo-{1,1}
    var bounds_hi = p_boxes_leaf[color].bounds.hi+{1,1}
    if bounds_lo.x < box_min[num_lvl] then
      bounds_lo.x = box_min[num_lvl]
    end
    if bounds_lo.y < box_min[num_lvl] then
      bounds_lo.y = box_min[num_lvl]
    end
    if bounds_hi.x > box_max[num_lvl] then
      bounds_hi.x = box_max[num_lvl]
    end
    if bounds_hi.y > box_max[num_lvl] then
      bounds_hi.y = box_max[num_lvl]
    end
--    c.printf("bounds_hi = (%i,%i)\n", bounds_hi.x, bounds_hi.y)
    c.legion_domain_point_coloring_color_domain(c_neighb, color, rect2d{bounds_lo, bounds_hi})
  end
  var p_boxes_neighb = partition(aliased, r_boxes_leaf, c_neighb, ispace(int2d, factorize(config.parallelism)))
  c.legion_domain_point_coloring_destroy(c_neighb)
  var p_particles_leaf = preimage(r_particles, p_boxes_leaf, r_particles.boxes)
  var p_particles_neighb = preimage(r_particles, p_boxes_neighb, r_particles.boxes)

  --
  -- Upward pass
  --
  var ts_start_up = c.legion_get_current_time_in_micros()
  c.printf("Starting P2M\n")

  -- Particle to multipole
  var token : int32 = 0
  __demand(__parallel)
  for color in p_boxes_leaf.colors do
    P2M(p_boxes_leaf[color], p_particles_leaf[color], num_lvl, p)
  end
  c.printf("P2M complete!\n")

  --Multipole to multipole
  c.printf("Starting M2M\n")
  for lvl = num_lvl-1, 0, -1 do
    var r_l = p_boxes_lvl[lvl]
    var r_c = p_boxes_lvl[lvl+1]
    var p_l = partition(equal,r_l,ispace(int2d,factorize(config.parallelism)))
    var p_c = partition(equal,r_c,ispace(int2d,factorize(config.parallelism)))
    for color in p_l.colors do
      M2M(p_l[color], p_c[color], p)
    end
  end
  c.printf("M2M complete!\n")

  var ts_stop_up = c.legion_get_current_time_in_micros()
  c.printf("Upward pass took %.4f sec\n", (ts_stop_up - ts_start_up) * 1e-6)

  --
  -- Downward pass
  --
  var ts_start_down = c.legion_get_current_time_in_micros()

  -- Multipole to local
  c.printf("Starting M2L\n")
  __demand(__parallel)
  for color in p_boxes_self.colors do
    M2L(p_boxes_self[color], p_boxes_ilist[color], p)
  end
  c.printf("M2L complete!\n")

  -- Local to local
  c.printf("Starting L2L\n")
  for lvl = 1, num_lvl do
    var r_l = p_boxes_lvl[lvl]
    var r_c = p_boxes_lvl[lvl+1]
    var p_l = partition(equal,r_l,ispace(int2d,factorize(config.parallelism)))
    var p_c = partition(equal,r_c,ispace(int2d,factorize(config.parallelism)))
    for color in p_l.colors do
      L2L(p_l[color], p_c[color], p)
    end
  end
  c.printf("L2L complete!\n")

  -- Force evaluation
  c.printf("Evaluating forces\n")
  __demand(__parallel)
  for color in p_boxes_leaf.colors do
    eval_forces(p_boxes_leaf[color], p_boxes_neighb[color], p_particles_leaf[color], p_particles_neighb[color], num_lvl, box_per_dim, p)
  end
  c.printf("Force evaluation complete!\n")

  var ts_stop_down = c.legion_get_current_time_in_micros()
  c.printf("Downward pass took %.4f sec\n", (ts_stop_down - ts_start_down) * 1e-6)

  -- Dump output
  for color in p_particles_leaf.colors do
    token = dump_forces(p_boxes_leaf[color], p_particles_leaf[color], config.output, token)
  end
--  dump_coeffs(r_boxes, config.output, p, token)

end

regentlib.start(toplevel)
