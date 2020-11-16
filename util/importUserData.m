function importUserData()
% Function to import the user data files in the input folder.
% The files imported are:
%   - airTemps.csv is imported to temp_matrix
%   - snow.csv is imported to snow_matrix
%   - salinity.csv is imported to salinity_matrix
%   - oceanFlux.csv is imported to ocean_matrix
%
% Zero snow is not allowed and is instead replaced with an
% arbitrarily small number, 10^-5 m here. This could be adjusted if
% numerical problems arise.


global snow_matrix temp_matrix salinity_matrix ocean_matrix lwrad_matrix cloud_matrix wind_matrix location_matrix pressure_matrix humidity_matrix lwrad_boolean Q_l_boolean Ri_b_matrix;


temp_matrix = readmatrix('input/airTemps.csv');

try
    snow_matrix = readmatrix('input/snow.csv');
catch
    snow_matrix = double.empty;
    disp("No snow depths specified. Defaulting to 0 m.")
end
try
    salinity_matrix = readmatrix('input/salinity.csv');
catch
    salinity_matrix = double.empty;
    disp("No salinity profile provided. Defaulting to constant salinity of 5 ppt.")
end
try
    ocean_matrix = readmatrix('input/oceanFlux.csv');
catch
    ocean_matrix = double.empty;
    disp("No ocean flux specified. Defaulting to 0 W/m^2 ocean flux.")
end
try
    lwrad_matrix = readmatrix('input/longWaveRadiation.csv');
    lwrad_boolean = 1;
catch
    lwrad_matrix = double.empty;
    lwrad_boolean = 0;
    disp("No longwave radiation specified. No long wave radiation interactions will be considered.")
end
try
    cloud_matrix = readmatrix('input/cloudiness.csv');
catch
    cloud_matrix = double.empty;
    disp("No cloudiness specified. Defaulting to 0.")
end
try
    wind_matrix = readmatrix('input/windSpeed.csv');
catch
    wind_matrix = double.empty;
    disp("No wind speed specified. Defaulting to 5 m/s.")
end
try
    location_matrix = readmatrix('input/location.csv');
catch
    location_matrix = double.empty;
    disp("No location specified. Defaulting to no solar radiation.")
end
try
    humidity_matrix = readmatrix('input/humidity.csv');
    Q_l_boolean = 1;
catch
    humidity_matrix = double.empty;
    Q_l_boolean = 0;
    disp("No specific humidity specified. Turbulant latent heat exchanges at the surface are set to 0 W/m^2.")
end
try
    pressure_matrix = readmatrix('input/airPressure.csv');
catch
    pressure_matrix = double.empty;
    disp("No air pressure specified. Defaulting to 101325 Pa (1 atm).")
end
try
    Ri_b_matrix = readmatrix('input/bulk_richardson.csv');
catch
    Ri_b_matrix = double.empty;
    disp("No bulk Richardson numbers specified. Defaulting to 0.")
end

%Check the files provided have been inputted correctly
if (~isempty(snow_matrix))
    if(length(snow_matrix(1, :)) ~= 2 || length(snow_matrix(:, 1)) ~= length(snow_matrix(:, 2)))
        error('The provided file snow.csv does not have the required structure.')
    end
end

if(length(temp_matrix(1, :)) ~= 2 || length(temp_matrix(:, 1)) ~= length(temp_matrix(:, 2)))
    error('The provided file airTemps.csv does not have the required structure.')
end

if (~isempty(salinity_matrix))
    if(length(salinity_matrix(1, :)) ~= 2 || length(salinity_matrix(:, 1)) ~= length(salinity_matrix(:, 2)))
        error('The provided file salinity.csv does not have the required structure.')
    end
end

if (~isempty(ocean_matrix))
    if(length(ocean_matrix(1, :)) ~= 2 || length(ocean_matrix(:, 1)) ~= length(ocean_matrix(:, 2)))
        error('The provided file oceanFlux.csv does not have the required structure.')
    end
end

if (~isempty(lwrad_matrix))
    if(length(lwrad_matrix(1, :)) ~= 2 || length(lwrad_matrix(:, 1)) ~= length(lwrad_matrix(:, 2)))
        error('The provided file longWaveRadiation.csv does not have the required structure.')
    end
end

if (~isempty(cloud_matrix))
    if(length(cloud_matrix(1, :)) ~= 2 || length(cloud_matrix(:, 1)) ~= length(cloud_matrix(:, 2)))
        error('The provided file cloudiness.csv does not have the required structure.')
    end
end
if (~isempty(wind_matrix))
    if(length(wind_matrix(1, :)) ~= 2 || length(wind_matrix(:, 1)) ~= length(wind_matrix(:, 2)))
        error('The provided file windSpeed.csv does not have the required structure.')
    end
end
if (~isempty(location_matrix))
    if(length(location_matrix(1, :)) ~= 3 || length(location_matrix(:, 1)) ~= length(location_matrix(:, 2)) || length(location_matrix(:, 1)) ~= length(location_matrix(:, 3)))
        error('The provided file location.csv does not have the required structure.')
    end
end
if (~isempty(pressure_matrix))
    if(length(pressure_matrix(1, :)) ~= 2 || length(pressure_matrix(:, 1)) ~= length(pressure_matrix(:, 2)))
        error('The provided file airPressure.csv does not have the required structure.')
    end
end
if (~isempty(humidity_matrix))
    if(length(humidity_matrix(1, :)) ~= 2 || length(humidity_matrix(:, 1)) ~= length(humidity_matrix(:, 2)))
        error('The provided file humidity.csv does not have the required structure.')
    end
end
if (~isempty(Ri_b_matrix))
    if(length(Ri_b_matrix(1, :)) ~= 2 || length(Ri_b_matrix(:, 1)) ~= length(Ri_b_matrix(:, 2)))
        error('The provided file bulk_richardson.csv does not have the required structure.')
    end
end

%Replace 0 m of snow with 10^-5 m
if ~isempty(snow_matrix)
    snow_z_array = snow_matrix(:, 2);
    snow_z_array(snow_z_array == 0) = 10^-5;
    time_array = snow_matrix(:, 1);
    
    snow_matrix = [time_array snow_z_array];
end

end