fID = fopen('test5.dat','w');

num_part = 48
num_dim = num_part^(1/2)
epsilon = 1e-5

fprintf(fID, 'NUM_PART %i\nMODEL_TYPE GRAVITY\n', num_part);
fprintf(fID, 'EPSILON %e\n', epsilon);

for ii = 1:num_part      
    fprintf(fID,'PART %f %f %f %f %f\n', randn, rand, rand, 0.0, 0.0);
end

fclose(fID);