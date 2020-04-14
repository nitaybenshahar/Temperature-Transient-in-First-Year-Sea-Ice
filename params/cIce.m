function c_ice = cIce(temp, salinity)
% Heat capacity of ice
%
% temp:     ice temperature measured in degrees C.
% salinity: sea ice salinity measured in grams/gram.

%The units for salinity in this formula are parts per thousand (ppt) and
%so we introduce the factor of 1000. Also c was assumed to be in units J/(g K),
%introducing the another factor of 1000 to convert to J/(kg K)

c_ice = 1000*(2.113 + 0.0075.*temp -0.0034*(salinity*1000) + 0.00008.*temp.*(salinity*1000) + 18.04*(salinity*1000)./(temp.^2));

end