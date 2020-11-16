function p = getAirPressure(t)
% A function to return the user specified air pressure at a certain time
%
% t:    desired time (seconds)
% 
% Note: pressure_matrix and time_offset must be specified globally before this function
% can be called


global pressure_matrix time_offset;

if ~isempty(pressure_matrix)
    data_time = (pressure_matrix(:, 1)-pressure_matrix(1,1))*24*60*60;      %in seconds
    p = interp1(data_time, pressure_matrix(:, 2), t + time_offset, 'linear', 101325);
else
    p = ones(length(t), 1)*101325;
end