function ze = calcZenith(t)

global start_date

cur_time = start_date + seconds(t);

location = getLocation(t);

[~, el] = SolarAzEl(reshape(cur_time, [], 1), location(1), location(2), 0);

el = max(0, el);

ze = deg2rad(90-el);

end