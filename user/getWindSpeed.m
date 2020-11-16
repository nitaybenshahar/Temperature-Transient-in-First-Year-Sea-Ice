function ws = getWindSpeed(t)
% A function to return the user specified windspeed at a certain
% time. (default 5m/s)
%
% t:    desired time (seconds)
% 
% Note: wind_matrix and time_offset must be specified globally before this function
% can be called


global wind_matrix time_offset;
if (~isempty(wind_matrix))
    data_time = (wind_matrix(:, 1)-wind_matrix(1,1))*24*60*60;      %in seconds

    ws = interp1(data_time, wind_matrix(:, 2), t + time_offset, 'linear', 5);
else
    ws = zeros(length(t), 1)+5;
end