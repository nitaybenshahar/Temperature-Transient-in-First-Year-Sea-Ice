function ri_b = getRi_b(t)
% A function to return the user specified bulk richerdson numbers at a certain time
%
% t:    desired time (seconds)
% 
% Note: Ri_b_matrix and time_offset must be specified globally before this function
% can be called


global Ri_b_matrix time_offset;
if ~isempty(Ri_b_matrix)
    data_time = (Ri_b_matrix(:, 1)-Ri_b_matrix(1,1))*24*60*60;      %in seconds
    ri_b = interp1(data_time, Ri_b_matrix(:, 2), t + time_offset, 'linear', 0);
else
    ri_b = zeros(length(t), 1);
end