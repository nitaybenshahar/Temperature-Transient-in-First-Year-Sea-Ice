function rad_heating = solarHeating(t, H_snow, H_ice, z_snow, z_ice, T_s)
% A function to return the solar heating at a certain time
%
% t:        desired time (seconds)
% H_snow:   thickness of the snow (m)
% H_ice:    depth of the ice (m)
% z_snow:   depth points to calculate solar heating at (in the snow)
% z_ice:    depth points to calculate solar heating at (in the ice)

global location_matrix

if isempty(location_matrix)
    rad_heating.ice = zeros(length(z_ice), 1);
    rad_heating.snow = zeros(length(z_snow), 1);
else
    %constants used:
    % snow light absorption decay rate
    K_snow = 50;
    %upper ice light absorption decay rate
    K1_ice = 20;
    %change over region
    z_c = 0.2;
    %lower ice light absorption decay rate
    K2_ice = 0.8;
    
    c = getCloudiness(t);
    snow_bool = H_snow > 10^-4;
    
    a = albedo(H_ice, H_snow, T_s);
    Q_0 = transpose((1-a).*(1-0.52*c).*surfaceFlux(t, T_s));
    
    rad_heating.snow = Q_0.*exp(-(H_snow-z_snow).*K_snow);
    
    if snow_bool
        rad_heating.ice = Q_0.*exp(-H_snow.*K_snow)*(heaviside(z_c-z_ice).*exp(-z_ice.*K1_ice)...
                        + heaviside(z_ice-z_c).*exp(-z_c.*K1_ice).*exp(-(z_ice-z_c).*K2_ice));
    else
        rad_heating.ice = Q_0.*(heaviside(z_c-z_ice).*exp(-z_ice.*K1_ice) ...
                        + heaviside(z_ice-z_c).*exp(-z_c.*K1_ice).*exp(-(z_ice-z_c).*K2_ice));
    end
end
end