fID = fopen('test1.dat','w');

num_part = 4^2
num_dim = num_part^(1/2)
epsilon = 1e-5

fprintf(fID, 'NUM_PART %i\nMODEL_TYPE GRAVITY\n', num_part);
fprintf(fID, 'EPSILON %e\n', epsilon);

for i = 0:1:int8(num_dim)-1
    for j = 0:1:int8(num_dim)-1
%             x = 1
%             y = 
%             z = 
        fprintf(fID,'PART %f %f %f %f %f\n', 10.0, i, j, 0.0, 0.0);
    end
end

fclose(fID);