function ocean_flux = getOceanFlux(t)
% A function to return the user specified oceanic heat flux at a given time
%
%   - t:    desired time (seconds)
%
% Note: ocean_matrix and time_offset must be specified before this function is called

global ocean_matrix time_offset
if (~isempty(ocean_matrix))
    times = ocean_matrix(:, 1)*24*60*60;    %in seconds
    fluxes = ocean_matrix(:, 2);
    ocean_flux = interp1(times, fluxes, t+time_offset, 'linear', 0);
else
    ocean_flux = zeros(length(t), 1);  %default to no ocean flux
end
end