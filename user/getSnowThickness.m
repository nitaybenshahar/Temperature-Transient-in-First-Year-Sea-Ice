function snow_thickness = getSnowThickness(t)
% A function to return the user specified snow thickness at a specific
% time.
%
% t:    desired time (seconds)
%
% Note: snow_matrix and time_offset must be specified globally before this
% function can be called
% If no snow has been inputted

global time_offset snow_matrix

if (~isempty(snow_matrix))
    snow_z_array = snow_matrix(:, 2);
    time_array = (snow_matrix(:, 1)-snow_matrix(1, 1))*60*60*24;    %to convert to seconds
    snow_thickness = interp1(time_array, snow_z_array, time_offset + t, 'linear', 0);
else
    snow_thickness = zeros(length(t), 1) + 10^-4;
end

end