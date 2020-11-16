function surface_temp = surfaceTemp(snow_temp, ice_temp, snow_thickness, ice_depth, T_s)
% A function to calculate the ice surface temperature between the interface of the snow
% and the ice as detailed in equations 10 of the write up.
%
% snow_temp:        temperature of the snow at the point of calculation closest
%                   to the snow/ice interface (ie at
%                   z=snow_thickness*delta_h_snow).
% ice_temp:         temperature of the ice at the point of calculation closest
%                   to the snow/ice interface (ie at
%                   z=-ice_thickness*delta_h_ice).
% snow_thickness:   thickness of the snow.
% ice_depth:        thickness of the ice.
% T_s:              upper surface temperature. This value is adopted if no
%                   snow is present.
%
%Note: delta_h_ice and delta_h_snow must be specified
%globally before this function is called.

global delta_h_ice delta_h_snow;

salinity = getSalinity(0);

k_ice = kIce(ice_temp, salinity);
k_snow = kSnow(snow_temp);


const_ice = -k_ice./(delta_h_ice*ice_depth);
const_snow = k_snow./(delta_h_snow*snow_thickness);

f_surface_temp = (const_ice.*ice_temp + const_snow.*snow_temp)./(const_ice + const_snow);

surface_temp(snow_thickness > 10^-4) = f_surface_temp(snow_thickness > 10^-4);
surface_temp(snow_thickness <= 10^-4) = ones(sum(snow_thickness <= 10^-4), 1).*T_s(snow_thickness <= 10^-4);


end