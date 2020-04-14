function deriv = iterateDEs(t, var)
    % Function to iterate the set of PDE's given in equations 23 and 30 in
    % the report.
    %
    % t:    desired time (seconds)
    % var:  a set of the temperature profiles (snow then ice) at the discretized positions along with the
    %       ice depth var=[snow_temp_profile ice_temp_profile ice_depth]
    
    %In this function, h_ice and h_snow is the reparamterized depth variable, living in
    %[0, 1]
    
    
    global num_of_points_ice num_of_points_snow T_freezing delta_h_ice delta_h_snow V_a;
    
    temp_snow = var(1:num_of_points_snow);
    temp_ice = var(num_of_points_snow+1:num_of_points_snow+num_of_points_ice);
    ice_depth = var(num_of_points_snow+num_of_points_ice+1);
    
    %fetch the air temperature
    air_temp = calcAirTemp(t);
    
    %snow thickness from the input data (user specified)
    snow_thickness = snowThickness(t);
    snow_growth_rate = snowGrowthRate(t);
    
    %initialize response variables
    dTice_dt = zeros(length(temp_ice),1);
    dTsnow_dt = zeros(length(temp_snow), 1);
    
    %Calculate the advancement of the boundary as per the boundary growth
    %equation
    %dT/dz=dT/dh dh/dz
    interface_salinity = calcSalinity(ice_depth);
    dT_dz = (T_freezing-temp_ice(end))/delta_h_ice * 1/ice_depth;
    ocean_flux = calcOceanFlux(t);
    growth_rate = iceGrowthConstant(T_freezing, interface_salinity, V_a)*dT_dz + ocean_flux/(LIce(T_freezing, interface_salinity) * rhoIce(T_freezing, interface_salinity, V_a));
    
    
    if snow_thickness > 0 %if snow is present
        T_snow_ice = surfaceTemp(temp_snow(1), temp_ice(1), snow_thickness, ice_depth, t);
    else
        T_snow_ice = air_temp;
    end
    
    %iterate the snow then ice temperature profile
    for i=1:length(temp_snow)
        h_snow = i/(length(temp_snow)+1);
        %equation 30
        if i==1
            dTsnow_dt(i) = alphaSnow(temp_snow(i))/(snow_thickness^2)*1/(delta_h_snow^2)*(temp_snow(i+1)-2*temp_snow(i)+T_snow_ice) + h_snow/snow_thickness *(temp_snow(i+1)-temp_snow(i))/delta_h_snow*snow_growth_rate;
        elseif i==length(temp_snow)
            dTsnow_dt(i) = alphaSnow(temp_snow(i))/(snow_thickness^2)*1/(delta_h_snow^2)*(air_temp-2*temp_snow(i)+temp_snow(i-1)) + h_snow/snow_thickness *(temp_snow(i)-temp_snow(i-1))/delta_h_snow*snow_growth_rate;
        else
            dTsnow_dt(i) = alphaSnow(temp_snow(i))/(snow_thickness^2)*1/(delta_h_snow^2)*(temp_snow(i-1)-2*temp_snow(i)+temp_snow(i+1)) + h_snow/snow_thickness *(temp_snow(i+1)-temp_snow(i-1))/(2*delta_h_snow)*snow_growth_rate;
        end
    end
    
    for i=1:length(temp_ice)
        h = i/(length(temp_ice)+1);
        %equation 30
        cp_z = cIce(temp_ice(i), calcSalinity(h*ice_depth)) * rhoIce(temp_ice(i), calcSalinity(h*ice_depth), V_a);
        if i==1
            k_ice_plus = kIce(0.5*(temp_ice(i)+temp_ice(i+1)), calcSalinity(ice_depth*(h+0.5)));
            k_ice_minus = kIce(0.5*(temp_ice(i)+T_snow_ice), calcSalinity(ice_depth*(h-0.5)));
            dTice_dt(i) = h/ice_depth *(temp_ice(i+1)-T_snow_ice)/(2*delta_h_ice)*growth_rate + 1/(cp_z)*1/(ice_depth^2)*1/(delta_h_ice^2)*(k_ice_minus*T_snow_ice-(k_ice_plus+k_ice_minus)*temp_ice(i)+k_ice_plus*temp_ice(i+1));
        elseif i==length(temp_ice)
            k_ice_plus = kIce(0.5*(temp_ice(i)+T_freezing), calcSalinity(ice_depth*(h+0.5)));
            k_ice_minus = kIce(0.5*(temp_ice(i)+temp_ice(i-1)), calcSalinity(ice_depth*(h-0.5)));
            dTice_dt(i) = h/ice_depth *(T_freezing-temp_ice(i-1))/(2*delta_h_ice)*growth_rate + 1/(cp_z)*1/(ice_depth^2)*1/(delta_h_ice^2)*(k_ice_minus*temp_ice(i-1)-(k_ice_plus+k_ice_minus)*temp_ice(i)+k_ice_plus*T_freezing);
        else
            k_ice_plus = kIce(0.5*(temp_ice(i)+temp_ice(i+1)), calcSalinity(ice_depth*(h+0.5)));
            k_ice_minus = kIce(0.5*(temp_ice(i)+temp_ice(i-1)), calcSalinity(ice_depth*(h-0.5)));
            dTice_dt(i) = h/ice_depth *(temp_ice(i+1)-temp_ice(i-1))/(2*delta_h_ice)*growth_rate + 1/(cp_z)*1/(ice_depth^2)*1/(delta_h_ice^2)*(k_ice_minus*temp_ice(i-1)-(k_ice_plus+k_ice_minus)*temp_ice(i)+k_ice_plus*temp_ice(i+1));
        end
    end
    
    %the response variable is the rate of change of each variable (the
    %snow/ice temperature profiles and the ice depth)
    deriv = [dTsnow_dt; dTice_dt; growth_rate];
end
