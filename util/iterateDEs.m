function deriv = iterateDEs(t, var)
    % Function to iterate the set of PDE's given in equations 7,8 and 9 in
    % the write up.
    %
    % t:    desired time (seconds)
    % var:  a set of the temperature profiles (snow then ice) at the discretized positions along with the
    %       ice depth. var=[snow_temp_profile ice_temp_profile ice_depth]
    
    %In this function, h_ice and h_snow is the reparamterized depth variable, living in
    %[0, 1]
    
    
    global num_of_points_ice num_of_points_snow T_freezing delta_h_ice delta_h_snow;
    
    temp_snow = var(1:num_of_points_snow);
    temp_ice = var(num_of_points_snow+1:num_of_points_snow+num_of_points_ice);
    ice_depth = var(num_of_points_snow+num_of_points_ice+1);
    
    %snow thickness from the input data (user specified)
    snow_thickness = getSnowThickness(t);
    snow_growth_rate = snowGrowthRate(t);
    
    
    %Calculate the advancement of the boundary as per the boundary growth
    %equation 9
    interface_salinity = getSalinity(ice_depth);
    dT_dz = (T_freezing-temp_ice(end))/delta_h_ice * 1/ice_depth;
    ocean_flux = getOceanFlux(t);
    growth_rate = kIce(T_freezing, interface_salinity)/(rhoIce(T_freezing, interface_salinity)...
        *LIce(T_freezing, interface_salinity))*dT_dz + ocean_flux/(LIce(T_freezing, interface_salinity) ...
        * rhoIce(T_freezing, interface_salinity));
    
    %If no snow is present ignore the effects of snow. Note the use of a
    %small parameter here as numerically we have trouble with exactly zero
    %snow.
    
    if snow_thickness > 10^-4
        T_1 = temp_snow(end);
        k_1 = kSnow(temp_snow(end));
        dz_1 = snow_thickness*delta_h_snow;
    else
        T_1 = temp_ice(1);
        k_1 = kIce(temp_ice(1), getSalinity(0));
        dz_1 = -ice_depth*delta_h_ice;
    end
    
    C_h = turb_heat(t);
    R = 8.3145;                 %universal gas constant
    M_a = M_air();              %molecular mass of air in kg/mol
    
    P = getAirPressure(t);
    H_air = getHumidity(t);
    T_air = getAirTemp(t);
    ws = getWindSpeed(t);
    
    rho_air = M_a.*P./(R.*(T_air+273.15).*(1+0.608.*H_air));
    T_s = upperTemp(k_1, dz_1, T_1, t, C_h.*rho_air, P, H_air, T_air, ws);
    T_snow_ice = surfaceTemp(temp_snow(1), temp_ice(1), snow_thickness, ice_depth, T_s);
    
    
    arr_h_snow = [0; transpose((1:length(temp_snow))*delta_h_snow); 1];
    arr_h_ice = [0; transpose((1:length(temp_ice))*delta_h_ice); 1];
    %calculate the solar heating at every required depth
    solar_heat = solarHeating(t, snow_thickness, -ice_depth, arr_h_snow*snow_thickness, -arr_h_ice*ice_depth, T_s);
    
    
    %iterate the snow then ice temperature profile according to equations 7
    %and 8
    %Note a slight deviation from the write up. The thermal conductivity of
    %sea ice/snow depends on the temperature and the salinity (i.e. the depth,
    %the salinity is calculated via the user specified salinity profile). 
    %Here to estimate the thermal conductivity in between sample points we 
    %linearly interpolate on the temperature and salinity and evaluate the 
    %thermal conductiviy at that point.
    
    %Calculate the change in the snow temperature profile
    %padded array values for calculation
    p_temp_snow = [T_snow_ice; temp_snow; T_s];
    p_temp_snow_s_minus = circshift(p_temp_snow, 1);
    p_temp_snow_s_plus = circshift(p_temp_snow, -1);
    
    arr_k_snow_plus = kSnow(0.5*(p_temp_snow+p_temp_snow_s_plus));
    arr_k_snow_minus = kSnow(0.5*(p_temp_snow+p_temp_snow_s_minus));
    
    arr_cp_snow = cSnow(p_temp_snow) .* rhoSnow();
    
    dTsnow_dt = 1./(arr_cp_snow*snow_thickness^2*delta_h_snow^2).*...
        (arr_k_snow_plus.*p_temp_snow_s_plus-(arr_k_snow_plus+arr_k_snow_minus)...
        .*p_temp_snow+arr_k_snow_minus.*p_temp_snow_s_minus) + arr_h_snow/snow_thickness...
        .*(p_temp_snow_s_plus-p_temp_snow_s_minus)/(2*delta_h_snow)*snow_growth_rate + solar_heat.snow./arr_cp_snow;
    %chop erroneous ends
    dTsnow_dt = dTsnow_dt(2:end-1);
    
    
    %Now apply a similar process to the ice temperatures
    
    %padded array values for calculation
    p_temp_ice = [T_snow_ice; temp_ice; T_freezing];
    p_temp_ice_s_minus = circshift(p_temp_ice, 1);
    p_temp_ice_s_plus = circshift(p_temp_ice, -1);
    
    arr_k_ice_plus = kIce(0.5*(p_temp_ice+p_temp_ice_s_plus), getSalinity((arr_h_ice+0.5*delta_h_ice)*ice_depth));
    arr_k_ice_minus = kIce(0.5*(p_temp_ice+p_temp_ice_s_minus), getSalinity((arr_h_ice-0.5*delta_h_ice)*ice_depth));
    
    arr_cp_ice = cIce(p_temp_ice, getSalinity(arr_h_ice*ice_depth)) .* rhoIce(p_temp_ice, getSalinity(arr_h_ice*ice_depth));
    dTice_dt = 1./(arr_cp_ice.*ice_depth^2.*delta_h_ice^2)...
        .*(arr_k_ice_minus.*p_temp_ice_s_minus-(arr_k_ice_plus+arr_k_ice_minus)...
        .*p_temp_ice+arr_k_ice_plus.*p_temp_ice_s_plus) + arr_h_ice.*growth_rate./ice_depth ...
        .*(p_temp_ice_s_plus-p_temp_ice_s_minus)/(2*delta_h_ice) + solar_heat.ice./arr_cp_ice;
    %chop erroneous ends
    dTice_dt = dTice_dt(2:end-1);
    
    deriv = [dTsnow_dt; dTice_dt; growth_rate];
end
