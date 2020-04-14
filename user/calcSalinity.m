function salinity = calcSalinity(depth)
% A function to return the user specified ice salinity (measured in grams/gram) at a certain ice depth.
% If no salinities are provided or the depth is outside the user specified
% range than the default salinity is 5/1000.
%
% depth:    specified ice depth in meters.
%
% Note: salinity_matrix must be specified before this function is called

global salinity_matrix
if (~isempty(salinity_matrix))
    depths = salinity_matrix(:, 1);
    salinities = salinity_matrix(:, 2);
    
    if (depth >= min(depths) && depth <= max(depths))
        salinity = interp1(depths, salinities, depth);
    else
        salinity = 5/1000;  %default salinity
    end
    
else
    salinity = 5/1000;  %default salinity
end

end