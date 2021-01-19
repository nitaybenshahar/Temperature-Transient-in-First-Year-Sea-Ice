function Q_s = surfaceFlux(t, T_s)
% A function to return the surface flux of solar radiation based on the
% solar azimuth.
%
% t:    desired time (seconds)

az = transpose(calcZenith(t));

%using Shine 1984 parameterization
S = 1362;           %solar constant
P = getAirPressure(t);
e = 6.1115*(1.0003+4.18*10^(-6)*P).*exp(22.452*T_s./(272.55+T_s));

Q_s = S.*cos(az).^2./(1.2*cos(az)+10^-3*e*(2.7+cos(az)) + 0.1);

end