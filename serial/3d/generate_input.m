fID = fopen('test2.dat','w');

num_part = 8^2
num_dim = num_part^(1/3)
epsilon = 1e-5

fprintf(fID, 'NUM_PART %i\nMODEL_TYPE GRAVITY\n', num_part);
fprintf(fID, 'EPSILON %e\n', epsilon);

for i = 0:1:int8(num_dim)-1
    for j = 0:1:int8(num_dim)-1
        for k = 0:1:int8(num_dim)-1
%             x = 1
%             y = 
%             z = 
            fprintf(fID,'PART %f %f %f %f %f %f %f\n', 10.0, i, j, k,...
                0.0, 0.0, 0.0);
        end
    end
end

fclose(fID);