function C_h = turb_heat(t)
% Calculates the sensible heat transfer based on the following inputs:


k = 0.4;    %von karman constant

%parameters as in Gryanik et al 2020.
e_m = 3*10^4;
e_t = 1.5*10^5;
Pr_0 = 0.98;
Ri_b = getRi_b(t);

Ri_b_bar = Ri_b/Pr_0;

A = (((log(e_m)+23.50).^5.25)./(181.3.*(log(e_t)+16.67).^2.625))...
    .*((((log(e_m)+23.50).^2)./(log(e_t)+16.67))-(log(e_m).^2)./(log(e_t)));

C_d = k.^2./((log(e_m)+50.0.*((1+0.3.*(log(e_m).^2./(log(e_t)).*Ri_b_bar +...
    A.*Ri_b_bar.^3.625)).^(1/3)-1)));
C_h = (k.*(C_d).^0.5)./(log(e_t)+12.5*log(1+0.4.*(log(e_m).^2./(log(e_t))...
    .*Ri_b_bar+A.*Ri_b_bar.^3.625)));

end