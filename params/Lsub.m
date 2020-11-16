function L = Lsub(temp)
% Latent heat of sublimation from ice/snow to air
%
% temp:     ice temperature measured in degrees C.

L = 10^5*(28.34-0.00149*(temp));

end