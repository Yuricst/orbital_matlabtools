function outputarr = cell2array(inputcell)
% function converts vertical cell to vertical array

outputarr = zeros(length(inputcell), 1);

for i = 1:length(inputcell)
    outputarr(i,1) = inputcell{i,1};
end

end