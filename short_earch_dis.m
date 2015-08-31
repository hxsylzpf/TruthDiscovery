function d = short_earch_dis(lat1, lon1, lat2, lon2)
% the input are decimal degrees
% need to convert to radians
lat1 = lat1/180*pi;
lon1 = lon1/180*pi;
lat2 = lat2/180*pi;
lon2 = lon2/180*pi;
R = 6371*10^3;

x = (lon2-lon1) * cos((lat1+lat2)/2);
y = (lat2-lat1);
d = sqrt(x*x + y*y) * R;
