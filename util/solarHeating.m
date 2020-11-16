function rad_heating = solarHeating(t, H_snow, H_ice, z_snow, z_ice, T_s)
% A function to return the solar heating at a certain time
%
% t:        desired time (seconds)
% H_snow:   depth of the snow (m)
% z_snow:   depth points to calculate solar heating at (in the snow)
% z_ice:    depth points to calculate solar heating at (in the ice)

global location_matrix

if isempty(location_matrix)
    rad_heating.ice = zeros(length(z_ice), 1);
    rad_heating.snow = zeros(length(z_snow), 1);
else
    %constants used:
    % snow light absorption decay rate
    v_snow = 1/100;
    %upper ice light absorption decay rate
    v_1_ice = 1/20;
    %change over region
    z_c = 0.2;
    %lower ice light absorption decay rate
    v_2_ice = 1/2;
    
    c = getCloudiness(t);
    snow_bool = H_snow > 10^-4;
    
    a = albedo(H_ice, H_snow, T_s);
    Q_0 = transpose((1-a).*(1-0.52*c).*surfaceFlux(t, T_s));
    
    rad_heating.snow = 1/v_snow*Q_0.*exp(-(H_snow-z_snow)/v_snow);
    
    if snow_bool
        rad_heating.ice = Q_0.*exp(-H_snow/v_snow)*(heaviside(z_c-z_ice).*exp(-z_ice/v_1_ice)/v_1_ice ...
                        + heaviside(z_ice-z_c).*exp(-z_c/v_1_ice).*exp(-(z_ice-z_c)/v_2_ice)/v_2_ice);
    else
        rad_heating.ice = Q_0.*(heaviside(z_c-z_ice).*exp(-z_ice/v_1_ice)/v_1_ice ...
                        + heaviside(z_ice-z_c).*exp(-z_c/v_1_ice).*exp(-(z_ice-z_c)/v_2_ice)/v_2_ice);
    end
end
end