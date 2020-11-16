function Q_lw = getLWrad(t)
% A function to return the user specified long wave radiation at a certain
% time. (default 0)
%
% t:    desired time (seconds)
% 
% Note: lwrad_matrix and time_offset must be specified globally before this function
% can be called


global lwrad_matrix time_offset;
if (~isempty(lwrad_matrix))
    data_time = (lwrad_matrix(:, 1)-lwrad_matrix(1,1))*24*60*60;      %in seconds

    Q_lw = interp1(data_time, lwrad_matrix(:, 2), t + time_offset, 'linear', 0);
else
    Q_lw = zeros(length(t), 1);
end