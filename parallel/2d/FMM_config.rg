import "regent"

local c = regentlib.c

struct FMMConfig
{
  input         : int8[512];
  output        : int8[512];
  num_particles : uint64;
  model_type    : int8[512];
  epsilon       : double;
  parallelism   : uint64;
}

local cstring = terralib.includec("string.h")

terra print_usage_and_abort()
  c.printf("Usage: regent.py FMM.rg [OPTIONS]\n")
  c.printf("OPTIONS\n")
  c.printf("  -h            : Print the usage and exit.\n")
  c.printf("  -i {file}     : Use {file} as input.\n")
  c.printf("  -o {file}     : Save the results to {file}.\n")
  c.abort()
end

terra file_exists(filename : rawstring)
  var file = c.fopen(filename, "rb")
  if file == nil then return false end
  c.fclose(file)
  return true
end

terra FMMConfig:initialize_from_command()
  var input_given = false
  var output_given = false
  var dummy_var : int8[512]
  self.epsilon = 1e-5
  self.parallelism = 1

  var args = c.legion_runtime_get_input_args()
  var i = 1
  while i < args.argc do
    -- Help
    if cstring.strcmp(args.argv[i], "-h") == 0 then
      print_usage_and_abort()

    -- Input file
    elseif cstring.strcmp(args.argv[i], "-i") == 0 then
      i = i + 1
      var file = c.fopen(args.argv[i], "rb")
      if file == nil then
        c.printf("File '%s' doesn't exist!\n", args.argv[i])
        c.abort()
      end
      cstring.strcpy(self.input, args.argv[i])
      c.fscanf(file, "%s %llu\n %s %s\n %s %lf",
               &dummy_var, &self.num_particles, &dummy_var, &self.model_type,
               &dummy_var, &self.epsilon)
      input_given = true
      c.fclose(file)

    -- Output file
    elseif cstring.strcmp(args.argv[i], "-o") == 0 then
      i = i + 1
      if file_exists(args.argv[i]) then
        c.printf("File '%s' already exists!\n", args.argv[i])
        c.abort()
      end
      cstring.strcpy(self.output, args.argv[i])
      output_given = true

    -- Help
    elseif cstring.strcmp(args.argv[i], "-p") == 0 then
      i = i + 1
      self.parallelism = c.atoi(args.argv[i])
    end

    i = i + 1
  end
  if not input_given then
    c.printf("Input file must be given!\n\n")
    print_usage_and_abort()
  end
  if not output_given then
    cstring.strcpy(self.output, "forces.out")
  end
end

return FMMConfig
