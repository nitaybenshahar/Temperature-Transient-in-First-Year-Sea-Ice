function model_output = seaIceModel()
% 1D sea ice thermodynamics model solved using the method of lines
% Snow depth is fixed externally.
%
% The model requires the following input:
%   - input/airTemps.csv with the atmospheric temperature as a function of
%   time.
%   - input/snow.csv with the snow thickness as a function of time.
%   - The simulation variables below.
%   - Optional: input/salinity.csv if a salinity profile is known.
%   - Optional: input/oceanFlux.csv if an oceanic flux is to be included.
%
%
% OUTPUT:
%   temp_profile_snow/ice:
%       - These are 2D arrays containing the calculated temperature
%       profiles of the snow/ice over time at the specified resolutions.
%
%   present_snow_thickness/present_ice_depth:
%       - These are 1D arrays containing the calculated snow/ice
%       thicknesses at every point in time. The snow thickness will simply
%       be an interpolation of the user defined snow thicknesses at the
%       desired times.
%
%   present_time/depth:
%       - 1D arrays containing the corresponding times/depth of the above
%       calculations
%
%   present_z_ice/snow:
%       - 1D arrays witht he corresponding depths/heights of the ice/snow
%       temperature profiles


%Get access to the necessary subroutines
addpath('params', 'user', 'util', 'input')

global num_of_points_ice num_of_points_snow delta_h_ice delta_h_snow time_offset T_freezing V_a snow_matrix temp_matrix;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%Simulation Variables to be adjusted%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%set up constants
num_of_points_ice = 30;
num_of_points_snow = 3;
odeset('RelTol',1e-9, 'AbsTol', 1e-10);                                      % The error tolerance of the ODE solver
final_t = 40.5*24*60*60;                                      % seconds
time_offset = 16.5*24*60*60;                                  % seconds
initial_ice_depth = -0.15;                                  % meters

%The output of the model will be in the following resolutions
%Note: This does not change the accuracy of the model. To get more accurate
%models increase num_of_points_ice/snow and/or decrease the error tolerance of the
%solver (ie odeset('RelTol', 1e-11); above)
present_dt = 60*60;
z_res_ice = 0.01;                                           % meters
z_res_snow = 0.01;                                          % meters


%Parameters
T_freezing = -1.87;                                         % Celsius
V_a = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read in the user provided air temperature, snow depth and salinity
% profile data
importUserData();

%Check the parameters allow the simulation to run
if time_offset + final_t > (temp_matrix(end, 1) - temp_matrix(1, 1))*24*60*60
    error('Air temperature data does not cover the requested simulation time.')
end
if time_offset + final_t > (snow_matrix(end, 1) - snow_matrix(1, 1))*24*60*60
    error('Snow cover data does not cover the requested simulation time.')
end
if initial_ice_depth >= 0
    error('The initial ice thickness must be negative (below sea level).')
end
if(num_of_points_ice < 1 || num_of_points_snow < 1)
    error('Not enough calculation points to run the simulation. Must have at least one point for both ice and snow, even if no snow is present')
end

delta_h_ice = 1/(num_of_points_ice+1);
delta_h_snow = 1/(num_of_points_snow+1);

%set up initial conditions
T_air = calcAirTemp(0);

initial_var = ones(num_of_points_ice+num_of_points_snow+1, 1)*T_air;
%[[Tsnow], [Tice], ice_depth]

%Initial temperature profile is set to be the steady state solution to the
%heat equation. Equation 24, 11 give us the ice/snow interface temperature.
%We assume a constant thermal conductivity for this estimation.
k_ice = kIce(T_freezing, calcSalinity(0));
k_snow = kSnow(T_freezing);
const_ice = -k_ice/initial_ice_depth;
const_snow = k_snow/snowThickness(time_offset);
T_ice_snow = (const_ice*T_freezing + const_snow*T_air)/(const_ice + const_snow);

for i = 1:num_of_points_snow
    initial_var(num_of_points_snow-i+1) = T_air + (T_ice_snow - T_air)/(num_of_points_snow+1)*i;
end
for i=1:num_of_points_ice
    initial_var(num_of_points_snow+i) = T_ice_snow + (T_freezing - T_ice_snow)/(num_of_points_ice+1)*i;
end

initial_var(end) = initial_ice_depth;
tspan = [0 final_t];

%Run ODE solver
[output_t, output_var] = ode15s(@(t, var) iterateDEs(t, var), tspan, initial_var, odeset);


%present the output as per the user selected data frequencies
present_time = ones(ceil(output_t(end)/present_dt), 1);
present_T_ice = ones(ceil(output_t(end)/present_dt), num_of_points_ice+2); %includes end points
present_T_snow = ones(ceil(output_t(end)/present_dt), num_of_points_snow+2); %includes end points
present_ice_depth = ones(ceil(output_t(end)/present_dt), 1);
present_snow_thickness = ones(ceil(output_t(end)/present_dt), 1);

for i=1:floor(output_t(end)/present_dt)
    time = (i-1)*present_dt;
    %interpolate on the variable list (temp profiles and ice depth)
    cur_var = interp1(output_t, output_var, time);
    my_T_snow = cur_var(1:num_of_points_snow);
    my_T_ice = cur_var(num_of_points_snow+1: num_of_points_snow+num_of_points_ice);
    my_ice_depth = cur_var(num_of_points_snow+num_of_points_ice+1);
    surf_temp = surfaceTemp(my_T_snow(1), my_T_ice(1), snowThickness(time), my_ice_depth, time);
    
    present_time(i) = time;
    present_T_ice(i, :) = [surf_temp my_T_ice T_freezing];
    present_T_snow(i, :) = [surf_temp my_T_snow calcAirTemp(time)];
    present_ice_depth(i) = my_ice_depth;
    present_snow_thickness(i) = snowThickness(time);
end


%present temperature profiles in real depth (non-frozen boundary)
temp_profile_ice = ones(length(present_time), ceil(max(-present_ice_depth)/z_res_ice))*T_freezing;
temp_profile_snow = zeros(length(present_time), ceil(max(present_snow_thickness)/z_res_snow));

present_z_ice = 0:-z_res_ice:min(present_ice_depth);
present_z_snow = 0:z_res_snow:max(present_snow_thickness);

for i=1:length(present_time)
    depths = 0:-z_res_ice:present_ice_depth(i);
    depth_index = floor(-present_ice_depth(i)/z_res_ice)+1;
    temp_profile_ice(i, 1:depth_index) = interp1(0:delta_h_ice*present_ice_depth(i):present_ice_depth(i), present_T_ice(i, :), depths);
    
    depths = 0:z_res_snow:present_snow_thickness(i);
    depth_index = floor(present_snow_thickness(i)/z_res_snow)+1;
    temp_profile_snow(i, 1:depth_index) = interp1(0:delta_h_snow*present_snow_thickness(i):present_snow_thickness(i), present_T_snow(i, :), depths);
end

model_output.time = present_time;
model_output.z_ice = present_z_ice;
model_output.z_snow = present_z_snow;
model_output.temp_profile_ice = temp_profile_ice;
model_output.temp_profile_snow = temp_profile_snow;
model_output.ice_depth = present_ice_depth;
model_output.snow_depth = present_snow_thickness;

end
