function surface_temp = surfaceTemp(snow_temp, ice_temp, snow_thickness, ice_depth, t)
% A function to calculate the ice surface temperature between the interface of the snow
% and the ice as detailed in equations 10 of the write up.
%
% snow_temp:        temperature of the snow at the point of calculation closest
%                   to the snow/ice interface (ie at z=snow_thickness*delta_h_snow)
% ice_temp:         temperature of the ice at the point of calculation closest
%                   to the snow/ice interface (ie at z=-ice_thickness*delta_h_ice)
% snow_thickness:   thickness of the snow
% ice_depth:        thickness of the ice
% t:                time (seconds), only used in the case that there is no
%                   snow for which we wish to set the surface to the air
%                   temperature
%
%Note: delta_h_ice and delta_h_snow must be specified
%globally before this function is called.

global delta_h_ice delta_h_snow;

salinity = calcSalinity(0);

k_ice = kIce(ice_temp, salinity);
k_snow = kSnow(snow_temp);

if snow_thickness > 0
    
    const_ice = -k_ice/(delta_h_ice*ice_depth);
    const_snow = k_snow/(delta_h_snow*snow_thickness);
    
    surface_temp = (const_ice*ice_temp + const_snow*snow_temp)/(const_ice + const_snow);
else
    surface_temp = calcAirTemp(t);
end

end