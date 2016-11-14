function M = GenerateRandomMeasurements( num_inputs, num_outputs )

% Generate random bases from the orthogonal colums of randomly 
% generated unitary matrices.         
M = {};
for j = 1:num_inputs
    U = RandomUnitary(num_outputs);
    for k = 1:num_outputs
        M{j}{k} = U(:,k)*U(:,k)';
		M{j}{k} = ( M{j}{k} + M{j}{k}' )/2;
    end
end    

end