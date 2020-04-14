function k_snow = kSnow(temp)
% Thermal conductivity of snow
%
% temp:     snow temperature measured in degrees C.


% The density in Yen 1981 was reported in Mg/m^3 so we divide the density by a factor
% of 1000 (converting from kg/m^3 to Mg/m^3).

k_snow = 0.0688 .* exp(0.0088.*temp + 4.6682.*rhoSnow()./1000);

end