function [new_filename]=ECEF2UTM(filename)
% ECEF to UTM
% Dette er en versjon av skriptet til Simen W, som er modifisert for aa 
% passe til import til andre program. All data presenteres i kolonner 
% med navn og i tillegg for koordinat har ogsaa standardavvik og kvalitet lagts til. 
% Kolonneformat: Punktnr., North, East, Kvalitet, Standardavvik NE
% Inn: .txt-fil ut: .csv-fil 
 
% Leser inn data ifra txt-fil til tabell
XYZ_format = readtable(filename);

% Quality og standard avvik beregnes direkt (std er i NE format (m))
Quality = XYZ_format.Var6; 
SDNE = sqrt(XYZ_format.Var8.^2 + XYZ_format.Var9.^2);

% Punktnummer
Point = [1:height(XYZ_format)]'; 

% Ellepsoide (wgs84)
a = 6378137.000000 ; b = 6356752.314245;

% Bruker ECEF2geode som tranformerer XYZ til Latitude longitude og høyde
Lat = zeros(height(XYZ_format)); Lon = zeros(height(XYZ_format)); 
for i = 1:height(XYZ_format)
    
    [Lat(i),Lon(i),h] = ECEF2geod(a, b, XYZ_format.Var3(i), XYZ_format.Var4(i), XYZ_format.Var5(i));
    
end

% Bruker deg2utm for å få Lat og long oversatt til UTM koordinater.
% Funksjonen er hentet ifra nettet.
[East, North, utmzone] = deg2utm(Lat(:,1), Lon(:,1));


% Setter in all data i en tabell
T2 = table(Point, North, East, Quality, SDNE);

% Skriver til csv format
split_filename = split(filename,"_");
split_filename = split(split_filename{3},".");
number = split_filename{1};
dato = strcat(number(7:8),number(5:6),number(1:4))
new_name = "CPOS_E_" + dato + ".csv";
writetable(T2, new_name,'Delimiter',',');

end

