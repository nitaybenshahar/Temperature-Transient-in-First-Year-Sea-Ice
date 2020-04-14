function K_ice = kIce(temp, salinity)
% Thermal conductivity of ice
%
% temp:     ice temperature measured in degrees C.
% salinity: sea ice salinity measured in grams/gram.
%
% Note: V_a must be specified globally before this function is called

global V_a

pure_ice_rho = 917;
sea_ice_rho = rhoIce(temp, salinity, V_a);

K_ice = sea_ice_rho./pure_ice_rho.*(2.11 - 0.011.*temp + 0.09*salinity./temp - (sea_ice_rho - pure_ice_rho)/1000); %Pringle 2006
end