function L_ice = LIce(temp, salinity)
% Heat of fusion of sea ice
%
% temp:     ice temperature measured in degrees C.
% salinity: sea ice salinity measured in grams/gram.

L_ice = (79.68 - 0.505.*temp - 27.3*salinity + 4311.5*salinity./temp)*4184;
%special case for no salinity where the freezing temp is 0C
L_ice(temp==0) = 79.68*4184;

end