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
Theta_T = @(T, i) T.*(P(i)./P_0).^(R./M_a./c_p_air());

h = C_h_rho_air.*ws.*c_p_air();

e = @(T, i) 6.1115*(1.0003 + 4.18*10^(-6)*P(i)).*exp(22.452*(T)./(272.55+T)); %T in Celsius
H_sat = @(T, i) 0.622*e(T, i)./P(i)./(1-0.378.*e(T, i)./P(i));                %T in Celsius
fun = @(T, i) (Q_d(i)-em*sigma.*T.^4)*lwrad_boolean ...
    + k(i).*(T_1(i)-T)./delta_z(i)+h(i).*(Theta_T(T_air, i)-Theta_T(T, i)) ...
    + Q_l_boolean*h(i).*Lsub(T-273.15)./c_p_air().*(H_air(i)-H_sat(T-273.15, i));

if length(T_1)>1
    %the following commented out code requires the Optimization Toolbox in Matlab. If
    %available, it may be used in replacement for the following for-loop
    %opts = optimoptions('fsolve', 'Display', 'off');
    %solns = fsolve(fun, T_1, opts);
    solns=273.15*ones(size(T_1));
    for i=1:length(T_1)
        solns(i) = fzero(@(T) fun(T, i), T_1(i));
    end
else
    solns = fzero(@(T) fun(T, 1), T_1);
end

T_s = solns-273.15;

end