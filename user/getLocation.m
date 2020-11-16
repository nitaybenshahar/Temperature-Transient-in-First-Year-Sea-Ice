function location = getLocation(t)
% A function to return the user specified location (latitude, longitude) at a certain time
%
% t:    desired time (seconds)
% 
% Note: temp_matrix and time_offset must be specified globally before this function
% can be called


global location_matrix time_offset;

if (~isempty(location_matrix))
    data_time = (location_matrix(:, 1)-location_matrix(1,1))*24*60*60;      %in seconds

    location = [interp1(data_time, location_matrix(:, 2), t + time_offset) ...
                interp1(data_time, location_matrix(:, 3), t + time_offset)];
else
    location = [nan nan];
end
end