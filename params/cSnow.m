function c_snow = cSnow(temp)
% Heat capacity of snow
%
% temp:     snow temperature measured in degrees C.

c_snow = (2.7442 + (temp+273.15).*0.1282)/18.02*1000;
end