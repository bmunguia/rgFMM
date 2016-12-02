--
-- Naive N-body simulation code which computes all interactions in O(N^2)
-- Authors : Brian C. Munguia
--

-- Import libraries
import "regent"
local c     = regentlib.c
local sqrt  = regentlib.sqrt(double)
local cmath = terralib.includec("math.h")
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
  field    : Vector2d; -- Field vector
  phi      : double;   -- Potential
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
  return c.fscanf(f, "%s %lf %lf %lf %lf %lf\n", &x, &source[0],
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

terra Vector2d.metamethods.__mul(v : Vector2d, c : double)
  return Vector2d { v._1 * c, v._2 * c }
end

terra comp_log(vec : Vector2d)
  return Vector2d {cmath.log(sqrt(cmath.pow(vec._1,2) + cmath.pow(vec._2,2))),
                   cmath.atan2(vec._2,vec._1)}
end

task initialize_particles(r_particles   : region(Particle),
                          filename      : int8[512])
where
  reads writes(r_particles)
do
  -- Set initial force and field to 0
  for particle in r_particles do
    particle.force = Vector2d {0.0, 0.0}
    particle.field = Vector2d {0.0, 0.0}
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

task dump_forces(r_particles : region(Particle),
                 filename : int8[512])
where
  reads(r_particles)
do
  var f = c.fopen(filename,"w")
  for particle in r_particles do
    c.fprintf(f, "%18.16f %18.16f %18.16f %18.16f %18.16f\n",
          particle.phi,
          particle.field._1, particle.field._2,
          particle.force._1, particle.force._2)
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

  -- Compute interactions
  var ts_start_forces = c.legion_get_current_time_in_micros()
  var i : uint8 = 0
  for particle1 in r_particles do
    for particle2 in r_particles do
      if particle1 ~= particle2 then
        var r : Vector2d = particle2.pos - particle1.pos
        var kern : double = r:kfun()
        var field_new : Vector2d = r * particle2.q * cmath.pow(kern,2)
        particle1.field = particle1.field + field_new
        particle1.force = particle1.force + field_new * particle1.q
	var phi_vec = comp_log(r) * particle2.q
	particle1.phi = particle1.phi + phi_vec._1
      end
    end
  end

  var ts_stop_forces = c.legion_get_current_time_in_micros()
  c.printf("Force calculation took %.4f sec\n", (ts_stop_forces - ts_start_forces) * 1e-6)

  dump_forces(r_particles, config.output)
end

regentlib.start(toplevel)
