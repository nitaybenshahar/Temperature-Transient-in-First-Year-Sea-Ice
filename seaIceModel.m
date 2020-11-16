function model_output = seaIceModel()
% 1D sea ice thermodynamics model solved for the provided forcing factors
% using the method of lines. See README file for required input and
% structure of the output.


%Get access to the necessary subroutines
addpath('params', 'user', 'util', 'input', 'SolarAzEl')

global num_of_points_ice num_of_points_snow delta_h_ice delta_h_snow time_offset T_freezing V_a temp_matrix start_date;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%Simulation Variables to be adjusted%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%set up constants
num_of_points_ice = 30;
num_of_points_snow = 20;
opts = odeset('RelTol',1e-6, 'AbsTol', 1e-4);                      % The error tolerance of the ODE solver
final_t = 80*24*60*60;                                      % seconds
time_offset = 0*24*60*60;                                   % seconds
initial_ice_depth = -0.15;                                  % meters

%The output of the model will be in the following resolutions
%This does not change the accuracy of the model. To get more accurate
%models increase num_of_points_ice/snow and/or decrease the error tolerance of the
%solver (ie odeset('RelTol', 1e-8); above)
present_dt = 60*60;
z_res_ice = 0.01;                                           % meters
z_res_snow = 0.001;                                          % meters

%for solar flux, need to define our time start
start_date = datetime('5/01/1997');


%Parameters
T_freezing = -1.8;                                         % Celsius
V_a = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read in the user provided air temperature, snow depth and salinity
% profile data
importUserData();

%Check the parameters allow the simulation to run
if time_offset + final_t > (temp_matrix(end, 1) - temp_matrix(1, 1))*24*60*60
    error('Air temperature data does not cover the requested simulation time.')
end
if initial_ice_depth >= 0
    error('The initial ice thickness must be negative (below sea level).')
end
if(num_of_points_ice < 1 || num_of_points_snow < 1)
    error('Not enough calculation points to run the simulation. Must have at least one point for both ice and snow, even if no snow is present')
end

%h represents the re-parameterized depth variable in the code
delta_h_ice = 1/(num_of_points_ice+1);
delta_h_snow = 1/(num_of_points_snow+1);

T_air = getAirTemp(0);

%set up initial conditions
initial_var = ones(num_of_points_ice+num_of_points_snow+1, 1)*T_air; %[[Tsnow], [Tice], ice_depth]

%Initial temperature profile is set to be the steady state solution to the
%heat equation assuming constant conductivity. Equation 10 gives us the ice/snow interface temperature.
%We assume a constant thermal conductivity for this estimation.
k_ice = kIce(T_freezing, getSalinity(0));
k_snow = kSnow(T_freezing);
const_ice = -k_ice/initial_ice_depth;
const_snow = k_snow/getSnowThickness(time_offset);
T_ice_snow = (const_ice*T_freezing + const_snow*T_air)/(const_ice + const_snow);


initial_var(1:num_of_points_snow) = T_air + (T_ice_snow - T_air)/(num_of_points_snow+1)*(num_of_points_snow:-1:1);
initial_var(num_of_points_snow+1:num_of_points_snow+num_of_points_ice) = T_ice_snow + (T_freezing - T_ice_snow)/(num_of_points_ice+1)*(1:num_of_points_ice);


initial_var(end) = initial_ice_depth;
tspan = 0:present_dt:final_t;

%Run ODE solver
[output_t, output_var] = ode15s(@(t, var) iterateDEs(t, var), tspan, initial_var, opts);

%the outputed solution is given in the user defined time spacing
present_time = output_t;

my_T_snow = output_var(:, 1:num_of_points_snow);
my_T_ice = output_var(:, 1+num_of_points_snow:num_of_points_snow+num_of_points_ice);
present_ice_depth = output_var(:, end);
present_snow_thickness = getSnowThickness(present_time);

%calculate upper layer temperatures
k=zeros(length(present_time), 1);
delta_z=zeros(length(present_time), 1);
T_1=zeros(length(present_time), 1);

k(present_snow_thickness>10^-4) = kSnow(my_T_snow(present_snow_thickness>10^-4, end));
k(present_snow_thickness<=10^-4) = kIce(my_T_ice(present_snow_thickness<=10^-4, 1), getSalinity(0));

delta_z(present_snow_thickness>10^-4) = delta_h_snow*present_snow_thickness(present_snow_thickness>10^-4);
delta_z(present_snow_thickness<=10^-4) = -delta_h_ice*present_ice_depth(present_snow_thickness<=10^-4);

T_1(present_snow_thickness>10^-4) = my_T_snow(present_snow_thickness>10^-4, end);
T_1(present_snow_thickness<=10^-4) = my_T_ice(present_snow_thickness<=10^-4, 1);

R = 8.3145; %universal gas constant
M_a = M_air();

P = getAirPressure(present_time);
H_air = getHumidity(present_time);
T_temp = getAirTemp(present_time);
ws = getWindSpeed(present_time);
rho_air = M_a.*P./(R.*(T_temp+273.15).*(1+0.608.*H_air));
C_h = turb_heat(present_time);

upper_temps = upperTemp(k, delta_z, T_1, present_time, C_h.*rho_air, P, H_air, T_air, ws);

%surf_temp refers to the surface between the snow and sea ice.
surf_temps = surfaceTemp(my_T_snow(:, 1), my_T_ice(:, 1), present_snow_thickness, present_ice_depth, upper_temps);

present_T_snow = [transpose(surf_temps) my_T_snow upper_temps];
present_T_ice = [transpose(surf_temps) my_T_ice ones(length(present_time), 1)*T_freezing];

present_z_ice = 0:-z_res_ice:min(present_ice_depth);
present_z_snow = 0:z_res_snow:max(present_snow_thickness);


real_depth_mesh_ice = present_ice_depth*(0:1/(num_of_points_ice+1):1);
real_depth_mesh_snow = present_snow_thickness*(0:1/(num_of_points_snow+1):1);


interp_ice = scatteredInterpolant(reshape(present_time*ones(1, length(0:1/(num_of_points_ice+1):1)), [],1), reshape(real_depth_mesh_ice, [],1), reshape(present_T_ice, [],1), 'linear', 'none');
[sampleX_i, sampleY_i] = meshgrid(present_z_ice, present_time);
temp_profile_ice = reshape(interp_ice(reshape(sampleY_i, [],1), reshape(sampleX_i, [],1)), length(present_time), length(present_z_ice));
temp_profile_ice(sampleX_i<present_ice_depth) = T_freezing;

interp_snow = scatteredInterpolant(reshape(present_time*ones(1, length(0:1/(num_of_points_snow+1):1)), [],1), reshape(real_depth_mesh_snow, [],1), reshape(present_T_snow, [],1), 'linear', 'none');
[sampleX_s, sampleY_s] = meshgrid(present_z_snow, present_time);
temp_profile_snow = reshape(interp_snow(reshape(sampleY_s, [],1), reshape(sampleX_s, [],1)), length(present_time), length(present_z_snow));
temp_profile_snow(sampleX_s>present_snow_thickness) = NaN;


%warn of times where melt would have occured
melt_index = sum([temp_profile_snow temp_profile_ice]>T_freezing, 2);

if sum(melt_index) ~= 0
    first_melt_index = find(melt_index, 1, 'first');
    warning("Melting temperatures occured at day "+num2str(present_time(first_melt_index)/60/60/24)+" of the simulation and possibly after.")
end

%Update response variables
model_output.time = present_time;
model_output.z_ice = present_z_ice;
model_output.z_snow = present_z_snow;
model_output.temp_profile_ice = temp_profile_ice;
model_output.temp_profile_snow = temp_profile_snow;
model_output.ice_depth = present_ice_depth;
model_output.snow_depth = present_snow_thickness;

end
