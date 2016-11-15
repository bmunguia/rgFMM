
--
-- The 'space_subdivision' task partitions the space into boxes.
-- For each box n we define the following subspaces:
--
task space_subdivision(r_particles : region(Particle))
where
  reads writes(r_particles)
do
  return 1
end

--
-- The 'R_expansion' task computes the potential phi_3^n(y) using the local
-- R expansion of phi(y,x_i) about the box center x_c^n.
-- The R expansion is given by:
--
--    phi_3^n(y) = sum_0^infty (A_m^n * R_m(y - x_c^n))
--      where A_m^n = sum_{x_i in I_3(n)} (q_i * a_m(x_i,x_c^n))
--            a_m(x_i,x_c) = -1/(x_i - x_c)^(m+1)
--            R_m(y - x_c) = (y - x_c)^m
--
task R_expansion(r_particles : region(Particle),
                 r_boxes     : region(Box(r_particles)))
where
  reads(r_particles,r_boxes), writes(r_boxes.phi)
do
  return 1
end

--
-- The 'S_expansion' task computes the potential phi_1^l(y) using the far field
-- S expansion of phi(y,x_i) about the box center x_c^n.
-- The S expansion is given by:
--
--    phi_1^l(y) = sum_0^infty (B_m^l * S_m(y - x_c^l))
--      where B_m^l = sum_{x_i in I_1(l)} (q_i * b_m(x_i,x_c^l))
--            b_m(x_i,x_c) = (x_i - x_c)^(m)
--            R_m(y - x_c) = 1/(y - x_c)^(m+1)
--
task S_expansion(r_particles : region(Particle),
                 r_boxes     : region(Box(r_particles)))
where
  reads(r_particles,r_boxes), writes(r_boxes.phi)
do
  return 1
end
