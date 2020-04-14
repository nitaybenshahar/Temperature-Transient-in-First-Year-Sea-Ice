function importUserData()
% Function to import the user data files in the input folder.
% The files imported are:
%   - airTemps.csv is imported to temp_matrix
%   - snow.csv is imported to snow_matrix
%   - salinity.csv is imported to salinity_matrix
%   - oceanFlux.csv is imported to ocean_matrix
%
% Zero snow is not allowed and is instead replaced with an
% arbitrarily small number, 10^-10m here. For a changing snow profile this
% should be set at 5*10^-4m to not cause trouble (found by trial and error)


global snow_matrix temp_matrix salinity_matrix ocean_matrix;

snow_matrix = readmatrix('input/snow.csv');
temp_matrix = readmatrix('input/airTemps.csv');
try
    salinity_matrix = readmatrix('input/salinity.csv');
catch
    salinity_matrix = double.empty;
    disp("No salinity profile provided. Defaulting to constant salinity of 5ppt")
end
try
    ocean_matrix = readmatrix('input/oceanFlux.csv');
catch
    ocean_matrix = double.empty;
    disp("No ocean flux specified. Defaulting to 0 ocean flux")
end

%Check the files provided have been inputted correctly
if(length(snow_matrix(1, :)) ~= 2 || length(snow_matrix(:, 1)) ~= length(snow_matrix(:, 2)))
    error('The provided file snow.csv does not have the required structure')
end

if(length(temp_matrix(1, :)) ~= 2 || length(temp_matrix(:, 1)) ~= length(temp_matrix(:, 2)))
    error('The provided file airTemps.csv does not have the required structure')
end

if (~isempty(salinity_matrix))
    if(length(salinity_matrix(1, :)) ~= 2 || length(salinity_matrix(:, 1)) ~= length(salinity_matrix(:, 2)))
        error('The provided file salinity.csv does not have the required structure')
    end
end

if (~isempty(ocean_matrix))
    if(length(ocean_matrix(1, :)) ~= 2 || length(ocean_matrix(:, 1)) ~= length(ocean_matrix(:, 2)))
        error('The provided file oceanFlux.csv does not have the required structure')
    end
end


%Replace 0m of snow with 10^-10m
snow_z_array = snow_matrix(:, 2);
snow_z_array(snow_z_array == 0) = 10^-10;
time_array = snow_matrix(:, 1);

snow_matrix = [time_array snow_z_array];

end