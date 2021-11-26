%Bruker deg2 utm og skriver ut i riktig format
filenames = ["20211014_123703_GMGD320.csv","20211019_125702.csv","20211021_110534.csv","20211026_124934_Siste runde.csv"];
for i = 1:length(filenames)
    filename = filenames(i)
    GPX_format_table = readtable(filename);
    GPX_format = [GPX_format_table.Var1, GPX_format_table.Var2, GPX_format_table.Var3];

    [x,y,utmzone] = deg2utm(GPX_format_table.Var2,GPX_format_table.Var1);

    split_filename = split(filename,"_");
    filename_dato = split_filename{1};
    dato = datetime(filename_dato,'InputFormat','yyyyMMdd');
    doy = day(dato,'dayofyear');
    current_year = year(dato);
    GPStime = current_year + doy/365;

    new_filename = filename_dato+".fri";
    
    utmzone

    fileID = fopen(new_filename,'wt');
    fprintf(fileID,'#PKT01 6615572.625   600117.213   136.039   2020/11/11\n');
    fprintf(fileID,'ITRF14-UTM(SONE-32),P,N,E,ELLITRF,EP,T\n');
    for i = 1:length(GPX_format)
        fprintf(fileID,'P%4.0f%14.3f%11.3f%10.3f%12.4f    \n',i+1000, y(i), x(i), GPX_format(i,3), GPStime); 
    end
    fclose(fileID);
end