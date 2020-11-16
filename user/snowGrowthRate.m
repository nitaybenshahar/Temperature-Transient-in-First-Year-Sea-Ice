function snow_growth_rate = snowGrowthRate(t)
% Returns the growth rate via interpolation of the snow depth via the user specified snow
% depths, at a specific time.
%
% t:    desired time (seconds)
%
% Note: snow_matrix and time_offset must be specified globally before this
% function can be called

global time_offset snow_matrix

if (~isempty(snow_matrix))
    snow_z_array = snow_matrix(:, 2);
    time_array = (snow_matrix(:, 1)-snow_matrix(1, 1))*60*60*24;    %to convert to seconds
    
    growth_rates = diff(snow_z_array)./(diff(time_array));
    growth_rates = [growth_rates; growth_rates(end)];
    snow_growth_rate = interp1(time_array, growth_rates, time_offset+t, 'previous', 0);
else
    snow_growth_rate = zeros(length(t), 1);
end

end