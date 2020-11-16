function T_s = upperTemp(k, delta_z, T_1, t, C_h_rho_air, P, H_air, T_air, ws)
%%Calculate the surface top surface temperature of snow/ice by balancing
%the incoming and outgoing heat at the interface.

global lwrad_boolean Q_l_boolean

%boltzman's constant
sigma = 1.38064852*10^(-23);

%emissivity
em = 0.99;

P_0 = 101325;               %reference pressure
R = 8.3145;                 %universal gas constant
M_a = M_air();              %molar mass of air, kg/mol


Q_d = getLWrad(t);
T_air = T_air+273.15;           %in K
T_1 = T_1 + 273.15;             %in K

%calculate potential temperatures
Theta_T = @(T) T.*(P./P_0).^(R./M_a./c_p_air());

h = C_h_rho_air.*ws.*c_p_air();

e = @(T) 6.1115*(1.0003 + 4.18*10^(-6)*P).*exp(22.452*(T)./(272.55+T)); %T in Celsius
H_sat = @(T) 0.622*e(T)./P./(1-0.378.*e(T)./P);                         %T in Celsius
fun = @(T) (Q_d-em*sigma*ones(length(k), 1).*T.^4)*lwrad_boolean ...
    + k.*(T_1-T)./delta_z+h.*(Theta_T(T_air)-Theta_T(T)) ...
    + Q_l_boolean*h.*Lsub(T-273.15)./c_p_air().*(H_air-H_sat(T-273.15));

if length(T_1)>1
    %the following requires the Optimization Toolbox in Matlab. Otherwise we
    %would have to apply fzero onto each row
    opts = optimoptions('fsolve', 'Display', 'off');
    solns = fsolve(fun, T_1, opts);
else
    solns = fzero(fun, T_1);
end

T_s = solns-273.15;

end