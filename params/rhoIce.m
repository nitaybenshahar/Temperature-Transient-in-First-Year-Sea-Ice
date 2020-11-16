function rho_ice = rhoIce(temp, salinity)
% Density of sea ice
%
% temp:     ice temperature measured in degrees C.
% salinity: sea ice salinity measured in grams/gram.
% 
% Note: V_a (the proportion of incorporated air in the ice, between 0 and
% 1) must be specified globally before this function is called

global V_a


%note the different units to Yen 1981, hence using 917 rather than 0.917

rho_ice = (1-V_a)*(1-4.51*salinity./temp)*917; 
end