function h = getHumidity(t)
% A function to return the user specified humidity at a certain time
%
% t:    desired time (seconds)
% 
% Note: humidity_matrix and time_offset must be specified globally before this function
% can be called


global humidity_matrix time_offset;
if ~isempty(humidity_matrix)
    data_time = (humidity_matrix(:, 1)-humidity_matrix(1,1))*24*60*60;      %in seconds
    h = interp1(data_time, humidity_matrix(:, 2), t + time_offset, 'linear', 10^-4);
else
    h = ones(length(t), 1)*10^-4;
end