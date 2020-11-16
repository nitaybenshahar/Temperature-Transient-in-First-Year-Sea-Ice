function salinity = getSalinity(z)
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
    salinity = interp1(depths, salinities, z, 'linear', 5/1000);
else
    salinity = ones(length(z), 1)*5/1000;  %default salinity
end

end