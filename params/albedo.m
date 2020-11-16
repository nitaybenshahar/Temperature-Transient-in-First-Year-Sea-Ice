function albedo = albedo(Hi, Hs, Ts)

global T_freezing

Hi = -Hi;

Tm = T_freezing;

%Here we use the CCSM3 parameterization of albedo
a_0 = 0.06;

a_i_vis = 0.73;
a_i_nir = 0.33;

a_s_vis = 0.96;
a_s_nir = 0.68;

a_it_vis = a_0 - a_i_vis.*min(atan(4*Hi)/atan(2), 1)+0.075.*min(Tm-Ts-1, 0);
a_it_nir = a_0 - a_i_nir.*min(atan(4*Hi)/atan(2), 1)+0.075.*min(Tm-Ts-1, 0);

a_vis = a_it_vis.*(1 - Hs/(Hs + 0.02)) + a_s_vis.*(Hs/(Hs + 0.02));
a_nir = a_it_nir.*(1 - Hs/(Hs + 0.02)) + a_s_nir.*(Hs/(Hs + 0.02));

albedo = 0.53*a_vis + 0.47*a_nir;

end