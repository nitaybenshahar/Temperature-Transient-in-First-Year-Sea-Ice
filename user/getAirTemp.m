function air_temp = getAirTemp(t)
% A function to return the user specified air temperature at a certain time
%
% t:    desired time (seconds)
% 
% Note: temp_matrix and time_offset must be specified globally before this function
% can be called


global temp_matrix time_offset;

data_time = (temp_matrix(:, 1)-temp_matrix(1,1))*24*60*60;      %in seconds

air_temp = interp1(data_time, temp_matrix(:, 2), t + time_offset);

end