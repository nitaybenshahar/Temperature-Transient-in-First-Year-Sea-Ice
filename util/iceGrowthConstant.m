function const = iceGrowthConstant(temp, salinity, V_a)
% Parameter to do with the moving sea ice boundary, as seen in equation 23 of the
% report
%
% temp:     Temperature at the ice/ocean interface
% salinity: Salinity at the ice/ocean interface
% V_a:      Proportion of incorporated air in the freezing sea ice

const = kIce(temp, salinity)/(rhoIce(temp, salinity, V_a)*LIce(temp, salinity));

end