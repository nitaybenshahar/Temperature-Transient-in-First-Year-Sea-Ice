function c = getCloudiness(t)
% A function to return the user specified cloudiness at a certain
% time. (default 0)
%
% t:    desired time (seconds)
% 
% Note: cloud_matrix and time_offset must be specified globally before this function
% can be called


global cloud_matrix time_offset;
if (~isempty(cloud_matrix))
    data_time = (cloud_matrix(:, 1)-cloud_matrix(1,1))*24*60*60;      %in seconds

    c = interp1(data_time, cloud_matrix(:, 2), t + time_offset, 'linear', 0);
else
    c = zeros(length(t), 1);
end