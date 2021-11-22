function [lat,lon,h]=ECEF2geod(a,b,X,Y,Z)

rho = sqrt((X^2)+(Y^2));

e2 = (a^2 - b^2)/a^2;

lat = atan(Z/rho);

d_lat = 1;


while d_lat > 1*10^(-11)
    
    N = a/(1-e2*sin(lat)^2)^(1/2);
    lat_new = atan((Z/rho)+(N*e2*sin(lat))/rho);
    d_lat = abs(lat_new - lat)*(180/pi);
    lat = lat_new;
    
end


lon = atan(Y/X);

lat = rad2deg(lat);
lon = rad2deg(lon);

h = (rho*cos(lat)) + (Z*sin(lat)) - (N*(1-e2*sin(lat)^2));

  
end