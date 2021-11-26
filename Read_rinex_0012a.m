% Isak / Tobias / Simen 
% Rinex Read
% Version 1.2
% Date: 09/11/2021

%%
% Store current directory
current_dir = pwd;

%readRinexObs304_dirpath = "C:\Users\tobia\OneDrive\Skrivebord\GMGD320\MATLAB\GNSS_reading_protocol\readRinexObs304\readRinexObs304_sourcecode";

%data dir relative to readRinexObs304_dirpath 
%data_dir   = "C:\Users\tobia\OneDrive\Skrivebord\GMGD320\New Folder\DATA\Topcon\Topcon";
%data_dir = 'C:\Users\isakf\Documents\Geomatikk\7_2021H\GMGD320\GMGD320_Group_Task\RINEX_Rounds';
Topcon_21   = '21_10_21_50002940.21o';
Topcon_19   = '19_10_21_50002920.21o';
Topcon_14   = '14_10_21_50002870.21o';
Topcon_26   = '26_10_21_50002990.21o';
Topcon_statisk = '320m3010 .21o';
Emlid_14 = 'reach_raw_202110141047.21O';
Emlid_19 = 'reach_raw_202110191050.21O';
Emlid_21 = 'reach_raw_202110210856.21O';
Emlid_statisk = 'reach_raw_202110280857.21O';
%filename    = 'reach_raw_202110210856.21O';
filename = append(pwd, '\RINEX_Rounds\', Topcon_statisk);
Topcon = "21o";
if contains(filename,Topcon) == 1
    system = "Topcon";
else 
    system = "Emlid" ;
end

%% read only GPS, with only code observation types and all bands. Don't read SS and LLI

% changing working directory to readRinexObs304
%cd(readRinexObs304_dirpath)

includeAllGNSSsystems = 0;
includeAllObsCodes = 0;
desiredGNSSsystems = ["G"];
desiredObsCodes = ["C", "L"];
desiredObsBands = [1,2,5];
readLLI  = 0;
readSS   = 0;

[GNSS_obs, GNSS_LLI, GNSS_SS, GNSS_SVs, time_epochs, nepochs, GNSSsystems,...
    obsCodes, approxPosition, max_sat, tInterval, markerName, rinexVersion, recType, timeSystem, leapSec, gnssType,...
    rinexProgr, rinexDate, antDelta, tFirstObs, tLastObs, clockOffsetsON, GLO_Slot2ChannelMap, success] = ...
    readRinexObs304(filename,readSS,readLLI,includeAllGNSSsystems,includeAllObsCodes, desiredGNSSsystems,... 
    desiredObsCodes, desiredObsBands);

%% Iterate through only GPS. For the first 10 epochs, compute the linear combination L1C-L2X if both L1C and L2X are present

if system == "Topcon" % Collects the correct phase and code code for Topcon
    desiredGNSS_systems = ["G"];
    phase1_code = "L1C";
    phase2_code = "L2W"; 
    code1_code  = "C1C";
    code2_code  = "C2W";
end

if system == "Emlid" % Collects the correct phase and code code for Emlid
    desiredGNSS_systems = ["G"];
    phase1_code = "L1C";
    phase2_code = "L2X";
    code1_code  = "C1C";
    code2_code  = "C2X";
end

c = 299792458; %Speed of light
phase1_wavelength = c/(1575.42*10^6);                   %Phase length for L1 in meter 
phase2_wavelength = c/(1227.60*10^6);                   %Phase length for L2 in meter 
alfa = (1575.42*10^6)^2/(1227.60*10^6)^2;               % L1/L2 in meters

n_epochs_to_display = nepochs;                          % Number of epochs

% map container to map from GNSS system code to actual name
GNSSnames = containers.Map(["G","R","C","E"], ["GPS", "GLONASS", "Beidou", "Galileo"]);
nGNSSsystems = length(desiredGNSS_systems);

MP1 = []; MP2 = [];                                     % Makes empty list to fill with the data of MP1 and MP2

for k = 1:nGNSSsystems
   %get current GNSS system
   sys = desiredGNSS_systems(k);
   
   % get current GNSS system index, return None if current sys is not in
   % data
   sysIndex = find([GNSSsystems{:}]==sys); 
   
   %fprintf('\nCurrently showing %s linear combinations\n\n', GNSSnames(sys))
   
   all_sat_MP1 = NaN(n_epochs_to_display,32);
   all_sat_MP2 = NaN(n_epochs_to_display,32);

   fasebrudd = zeros(n_epochs_to_display,32);
   IOD_indicator = NaN(n_epochs_to_display,32);
   % iterate through epochs for the current GNSS system
   for epoch =1:n_epochs_to_display
      % get number of sat for current GNSS system that have observations
      %Number of satelites for the corrent epoch.
      n_sat = GNSS_SVs{sysIndex}(epoch, 1);
      
      % get sat IDs of current GNSS system with observation this epoch 
      SVs = GNSS_SVs{sysIndex}(epoch, 2:n_sat+1);
      
      %iterate through all satelites of current GNSS system this epoch
      for sat = 1:n_sat
         SV = SVs(sat);
                
         % get index of phase 1 & 2 and code 1 & 2 obs for current sat if present
         phase1_index = ismember(obsCodes{sysIndex},phase1_code);
         phase2_index = ismember(obsCodes{sysIndex},phase2_code);
         code1_index  = ismember(obsCodes{sysIndex},code1_code);
         code2_index  = ismember(obsCodes{sysIndex},code2_code);
         
         % get phaseobservations and convert from units of cycles to meters
         phase1 = GNSS_obs{sysIndex}(SV, phase1_index, epoch)*phase1_wavelength;
         phase2 = GNSS_obs{sysIndex}(SV, phase2_index, epoch)*phase2_wavelength;
         
         % get code observations (m)
         code1 = GNSS_obs{sysIndex}(SV, code1_index, epoch);
         code2 = GNSS_obs{sysIndex}(SV, code2_index, epoch);
         if ~any([phase1, phase2] == 0)
            % Fasebruddindikator
            IOD = (alfa/(alfa-1))*(phase1 - phase2); 
            fasebrudd(epoch, SV) = IOD;
            IOD_indicator(epoch, SV) = IOD;
         end
         %check that none of the phase obs are missing
         if ~any([phase1, phase2, code1] == 0)
            % Multipath + Bias
            Mp1 = code1 - (1 + 2/(alfa-1))*phase1 + (2/(alfa-1))*phase2;                   
            all_sat_MP1(epoch,SV) = Mp1;
         end 
         if ~any([phase1, phase2, code2] == 0)
            Mp2 = code2 - (2*alfa/(alfa-1))*phase1 + (2*alfa/(alfa-1)-1)*phase2;                      
            all_sat_MP2(epoch,SV) = Mp2;
           
         end   
      end
   end
end
% Make the matrix that has only all the phaceshifts, and not other changes.
% so that you can plot the phaseshifts on its own
fasebrudd_plot = [];
[m,n] = size(IOD_indicator);
for i = 1:m-1
    for j = 1:n
        d2 = abs(IOD_indicator(i+1,j));
        d = abs(IOD_indicator(i,j));
        dis = d2 - d;
        if abs(dis)  > 4/60 ||  (isnan(d2) && ~isnan(d))
            if any(~isnan(all_sat_MP1(i+1:nepochs,j)))
                fasebrudd_plot = [fasebrudd_plot;j,i];
            end
        end 
    end
end 


% Make matrix that said if there is changes big enough to be considered as
% phaseshifts
[m,n] = size(fasebrudd);
d = [];
index = [];
for i = 1:m-1
    for j = 1:n
        d2 = fasebrudd(i+1,j);
        d = fasebrudd(i,j);
        dis = abs(d2) - abs(d);
        if abs(dis)  > 4/60 || (dis ~= 0 && i == 1)
            index = [index;j,i];
            %index = [index;j,i+1];
        elseif abs(dis)  < 4/60 && (dis ~= 0 && i+1 == nepochs)
            index = [index;j,i+1];
            %index = [index;j,i+1];
        end
    end
end  

index = sortrows(index,1);
[n, m] =size(index);
% Mean table is a table that has mean elements of the periods in it
%[satelite,first element in period, last element, mean value of the period]

mean_table_MP1 = [];
mean_table_MP2 = [];
for k = 1:n-1
    if (index(k,1) == index(k+1,1))% && index(k,2) ~= 1   %&& (k+2 < m && isnan(all_sat_MP1(index(k+2,1),index(k+2,2))) == 1)
        mean_table_MP1 = [mean_table_MP1;index(k,1),index(k,2),(index(k+1,2)),mean(all_sat_MP1(index(k,2)+1:(index(k+1,2)),index(k,1)),'omitnan')];
        mean_table_MP2 = [mean_table_MP2;index(k,1),index(k,2),(index(k+1,2)),mean(all_sat_MP2(index(k,2)+1:(index(k+1,2)),index(k,1)),'omitnan')];
    end
end

mean_table_MP1( any( isnan(mean_table_MP1), 2 ), : ) = [];
mean_table_MP2( any( isnan(mean_table_MP2), 2 ), : ) = [];
[m,n] = size(mean_table_MP1);
for k = 1:m-1
    if mean_table_MP1(k+1,2) == mean_table_MP1(k,3)
        mean_table_MP1(k+1,2) = mean_table_MP1(k+1,2)+1;
        %mean_table_MP1(k+1,4) = mean(all_sat_MP1(index(k+1,2)+1:(index(k+2,2)),index(k,1)),'omitnan');
    end
    if mean_table_MP2(k+1,2) == mean_table_MP2(k,3) 
        mean_table_MP2(k+1,2) = mean_table_MP2(k+1,2)+1;
    end
end

%% IOD 
[j,f] = size(mean_table_MP1);
for k = 1:j
    IOD_indicator(mean_table_MP1(k,2):mean_table_MP1(k,3),mean_table_MP1(k,1)) = IOD_indicator(mean_table_MP1(k,2):mean_table_MP1(k,3),mean_table_MP1(k,1))- IOD_indicator(mean_table_MP1(k,2),mean_table_MP1(k,1));
end
[n,k] = size(IOD_indicator);
IOD = figure
satelites = unique(mean_table_MP1(:,1));
legend_list = append('PRN',string(satelites(:)));
value = 0;
for i = 1:k
    plot(IOD_indicator(:,i))
    hold on
end
tekst = [system, 'Estimated ionospheric delay on the F1 frequensy'];
title(tekst);
ylim([-10 10]);
xlim([0 7700]);
xlabel('Epochs') 
x0=10;
y0=10;
width=800;
height=500;
set(gcf,'position',[x0,y0,width,height])
legend(legend_list,'Location','southeast','NumColumns',4)
ylabel('Noise(meters)')
filnavn = append('IOD_',system,'.png');
exportgraphics(IOD,filnavn)
%plot(real_MP(:,4))

%% Plot that shows MP1 and MP2 with just substraction of the first value.

% Fjerner Bias fra MP (MP-Bias)    
[m,n] = size(mean_table_MP1);
real_MP1 = all_sat_MP1;
real_MP2 = all_sat_MP2;
for i = 1:m
        real_MP1(mean_table_MP1(i,2):mean_table_MP1(i,3),mean_table_MP1(i,1)) = real_MP1(mean_table_MP1(i,2):mean_table_MP1(i,3),mean_table_MP1(i,1)) - real_MP1(mean_table_MP1(i,2)+1,mean_table_MP1(i,1));
        real_MP2(mean_table_MP2(i,2):mean_table_MP2(i,3),mean_table_MP2(i,1)) = real_MP2(mean_table_MP2(i,2):mean_table_MP2(i,3),mean_table_MP2(i,1)) - real_MP2(mean_table_MP2(i,2)+1,mean_table_MP2(i,1));
end
    
%Fjerner alle kolonner med kun "NaN" verdier til å plotte
plot_MP1 = real_MP1(:,~all(isnan(real_MP1)));
plot_MP2 = real_MP2(:,~all(isnan(real_MP2)));
satelites = unique(mean_table_MP1(:,1));
legend_list = append('PRN',string(satelites(:)));
% Plot
[m, n] = size(plot_MP1);
[m, n] = size(plot_MP2);

% Test

MP1 = figure;
for i = 1:n
    plot(plot_MP1(:, i))
    hold on
end
[row,kolonne] =  size(plot_MP1);
mean_value = mean(abs(plot_MP1(1:row,kolonne)),'omitnan');
tekst = [system, 'Estimated multipath on the P1-code(MP1) with substraction of the first element in the phase corrected arc of observations '];
title(tekst);
ylim([-10 10]);
xlim([-5 7700]);
x0=10;
y0=10;
width=1200;
height=500;
set(gcf,'position',[x0,y0,width,height])
legend(legend_list,'Location','southeast','NumColumns',4)
xlabel('Epochs') 
ylabel('Noise(meters)') 
filnavn = append('MP1_sub_',system,'.png');
exportgraphics(MP1,filnavn)
MP2 = figure;


for i = 1:n
    plot(plot_MP2(:, i))
    hold on
    
end
[row,kolonne] =  size(plot_MP2);
mean_value = mean(abs(plot_MP2(1:row,kolonne)),'omitnan');
tekst = [system, 'Estimated multipath on the P2-code(MP2) with substraction of the first element in the phase corrected arc of observations '];
title(tekst);
ylim([-10 10]);
xlim([-5 7700]);
xlabel('Epochs') 
x0=10;
y0=10;
width=1200;
height=500;
set(gcf,'position',[x0,y0,width,height])
legend(legend_list,'Location','southeast','NumColumns',4)
ylabel('Noise(meters)')
filnavn = append('MP2_sub_',system,'.png');
exportgraphics(MP2,filnavn)

%% Plot that shows MP1 and MP2 with just substraction of the mean value of the phase corrected arc of observations.


% Fjerner Bias fra MP (MP-Bias)    
[m,n] = size(mean_table_MP1);
real_MP1 = all_sat_MP1;
real_MP2 = all_sat_MP2;
for i = 1:m
    real_MP1(mean_table_MP1(i,2):mean_table_MP1(i,3),mean_table_MP1(i,1)) = (real_MP1(mean_table_MP1(i,2):mean_table_MP1(i,3),mean_table_MP1(i,1))) - (mean_table_MP1(i,4));
    real_MP2(mean_table_MP2(i,2):mean_table_MP2(i,3),mean_table_MP2(i,1)) = (real_MP2(mean_table_MP2(i,2):mean_table_MP2(i,3),mean_table_MP2(i,1))) - (mean_table_MP2(i,4));
end
    
%Fjerner alle kolonner med kun "NaN" verdier til å plotte
plot_MP1 = real_MP1(:,~all(isnan(real_MP1)));
plot_MP2 = real_MP2(:,~all(isnan(real_MP2)));
satelites = unique(mean_table_MP1(:,1));
legend_list = append('PRN',string(satelites(:)));
% Plot
[m, n] = size(plot_MP1);
[m, n] = size(plot_MP2);

% Test

MP1 = figure;
for i = 1:n
    plot(plot_MP1(:, i))
    hold on
end
[row,kolonne] =  size(plot_MP1);
mean_value = mean(abs(plot_MP1(1:row,kolonne)),'omitnan');
tekst = [system, 'Estimated multipath on the P1-code(MP1) with with the mean of the phase corrected arc of observations '];
title(tekst);
ylim([-10 10]);
xlim([-5 7700]);
x0=10;
y0=10;
width=1200;
height=500;
set(gcf,'position',[x0,y0,width,height])
legend(legend_list,'Location','southeast','NumColumns',4)
xlabel('Epochs') 
ylabel('Noise(meters)') 
filnavn = append('MP1_',system,'.png');
exportgraphics(MP1,filnavn)
MP2 = figure;
for i = 1:n
    plot(plot_MP2(:, i))
    hold on
end
[row,kolonne] =  size(plot_MP2);
mean_value = mean(abs(plot_MP2(1:row,kolonne)),'omitnan');
tekst = [system, 'Estimated multipath on the P2-code(MP2) with the mean of the phase corrected arc of observations '];
title(tekst);
x0=10;
y0=10;
width=1200;
height=500;
set(gcf,'position',[x0,y0,width,height])
ylim([-10 10]);
xlim([-5 7700]);
xlabel('Epochs') 
legend(legend_list,'Location','southeast','NumColumns',4)
ylabel('Noise(meters)')
filnavn = append('MP2_',system,'.png');
exportgraphics(MP2,filnavn)
%plot(real_MP(:,4))

[row,kolonne] =  size(mean_table_MP1);
[nr_sat,m] = size(satelites)
bias_table_MP1 = zeros(31,3);
bias_table_MP2 = zeros(31,3);

for k = 1:row
    bias_table_MP1(mean_table_MP1(k,1),:) = [mean_table_MP1(k,1),bias_table_MP1(mean_table_MP1(k,1),2)+ mean_table_MP1(k,3)-mean_table_MP1(k,2),bias_table_MP1(mean_table_MP1(k,1),3)+1];
    bias_table_MP2(mean_table_MP2(k,1),:) = [mean_table_MP2(k,1),bias_table_MP2(mean_table_MP2(k,1),2)+ mean_table_MP2(k,3)-mean_table_MP2(k,2),bias_table_MP2(mean_table_MP2(k,1),3)+1];
end
bias_table_MP2 = bias_table_MP2(any(bias_table_MP2,2),:);
bias_table_MP1 = bias_table_MP1(any(bias_table_MP1,2),:);

S_mp1 = [];
S_mp2 = [];
sum_mp = [];
m_mp1 = []; % mean mp1 
m_mp2 = []; % mean mp1 
for k = 1:nr_sat
    
    rsum_mp1 = sum(abs(real_MP1(:,satelites(k))),'omitnan'); % sum mp1
    rsum_mp2 = sum(abs(real_MP2(:,satelites(k))),'omitnan'); % sum mp2
    sum_mp1 = sum(real_MP1(:,satelites(k)).^2,'omitnan');
    sum_mp2 = sum(real_MP2(:,satelites(k)).^2,'omitnan');
    nevner_mp1 = (bias_table_MP1(k,2)-bias_table_MP1(k,3));
    nevner_mp2 = (bias_table_MP2(k,2)-bias_table_MP2(k,3));
    m_mp1 = [m_mp1;satelites(k),rsum_mp1/nevner_mp1]; % mean mp1
    m_mp2 = [m_mp2;satelites(k),rsum_mp2/nevner_mp2]; % mean mp2
    S_mp1 = [S_mp1;satelites(k),sqrt(sum_mp1/nevner_mp1)];
    S_mp2 = [S_mp2;satelites(k),sqrt(sum_mp2/(nevner_mp2))];
end 


%% Compute the RMS values for MP1 and MP2 

S_MP = figure;
S = [S_mp1(:,2),S_mp2(:,2)]; % Isak orginal
%S = [S_mp1(:,2),m_mp1(:,2)]; % Tobias kopi... 
h = bar(S_mp1(:,1),S,1);
hold on
mean_S_mp1 = mean(S_mp1(:,2)) ;
plot(xlim,[mean_S_mp1 mean_S_mp1])
hold on
mean_S_mp2 = mean(S_mp2(:,2)) ;
plot(xlim,[mean_S_mp2 mean_S_mp2])

tekst = ['RMS value for MP1 and MP2', system,]; 
legend_name = {'MP1','MP2','Mean MP1','Mean MP2'};
%tekst = ['Standard deviation for MP1 and MP2, for ', system,'mean MP1', mean(S_mp1(:,2)),'mean MP2', mean(S_mp2(:,2))];
% set 3 display names for the 3 handles
%set(h, {'DisplayName'}, {'MP1','MP2'}')
% Legend will show names for each color
x0=10;
y0=10;
width=1200;
height=500;
set(gcf,'position',[x0,y0,width,height])
legend(legend_name,'Location','northeast','NumColumns',2)

title(tekst);
xlabel('Satelite number') 
ylabel('standard deviation in meters')
filnavn = append('S_MP_',system,'.png');

exportgraphics(S_MP,filnavn)

%% Information on the detected phase-breach

% Number of sateliteobservations.
[n,l] = size(real_MP1);
i = 0;
for k = 1:n
    for j = 1:l
        if ~isnan(real_MP1(k,j))
            i = i+1;
        end
    end
end
[n,j] = size(fasebrudd_plot)
nr_fasebrudd = zeros(31,2);
for k = 1:n
    nr_fasebrudd(fasebrudd_plot(k,1),:) = [fasebrudd_plot(k,1), (nr_fasebrudd(fasebrudd_plot(k,1),2)+1)];
end
nr_fasebrudd = nr_fasebrudd(any(nr_fasebrudd,2),:); 

fasebrudd = figure;
h = bar(nr_fasebrudd(:,1),nr_fasebrudd(:,2),1);

tekst = "Number of phaseshift per satelite";

x0=10;
y0=10;
width=1200;
height=500;
set(gcf,'position',[x0,y0,width,height])
legend('Phaseshifts','Location','northwest','NumColumns',2)
title(tekst);
xlabel('Satelite number') 
ylabel('Nr of phaseshifts')
filnavn = append('phaseshifts_',system,'.png');
exportgraphics(fasebrudd,filnavn)
