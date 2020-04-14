function ocean_flux = calcOceanFlux(t)
% A function to return the user specified oceanic heat flux at a given time
%
%   - t:    desired time (seconds)
%
% Note: ocean_matrix and time_offset must be specified before this function is called

global ocean_matrix time_offset
if (~isempty(ocean_matrix))
    times = ocean_matrix(:, 1)*24*60*60;    %in seconds
    fluxes = ocean_matrix(:, 2);
    
    if (t+time_offset >= min(times) && t+time_offset <= max(times))
        ocean_flux = interp1(times, fluxes, t+time_offset);
    else
        ocean_flux = 0;  %default to no ocean flux
    end
    
else
    ocean_flux = 0;  %default to no ocean flux
end

end