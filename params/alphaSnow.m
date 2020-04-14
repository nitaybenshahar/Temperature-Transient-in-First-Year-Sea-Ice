function alpha_snow = alphaSnow(temp)
% Thermal diffusivity of snow is calculated via alpha=k/(rho*c)
%
% temp:     temperature of the snow in degrees C.

K = kSnow(temp);
c = cSnow(temp);
rho_snow = rhoSnow();

alpha_snow = K./(rho_snow.*c);
end