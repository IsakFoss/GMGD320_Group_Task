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
data_dir = "C:\Users\tobia\OneDrive\Skrivebord\GMGD320\New Folder\DATA\reach_raw_202110280857_RINEX_3_03 (1)";
%filename    = "320m3010.21o";
filename    = "reach_raw_202110280857.21O";
filename = append(data_dir, '\', filename);

%% read only GPS, with only code observation types and all bands. Don't read SS and LLI

% changing working directory to readRinexObs304
%cd(readRinexObs304_dirpath)

includeAllGNSSsystems = 0;
includeAllObsCodes = 0;
desiredGNSSsystems = ["G"];
desiredObsCodes = ["C", "L"];
desiredObsBands = [1,2];
readLLI  = 0;
readSS   = 0;

[GNSS_obs, GNSS_LLI, GNSS_SS, GNSS_SVs, time_epochs, nepochs, GNSSsystems,...
    obsCodes, approxPosition, max_sat, tInterval, markerName, rinexVersion, recType, timeSystem, leapSec, gnssType,...
    rinexProgr, rinexDate, antDelta, tFirstObs, tLastObs, clockOffsetsON, GLO_Slot2ChannelMap, success] = ...
    readRinexObs304(filename,readSS,readLLI,includeAllGNSSsystems,includeAllObsCodes, desiredGNSSsystems,... 
    desiredObsCodes, desiredObsBands);

%% Iterate through only GPS. For the first 10 epochs, compute the linear combination L1C-L2X if both L1C and L2X are present

desiredGNSS_systems = ["G"];
phase1_code = "L1C";
phase2_code = "L2X";
code1_code  = "C1C";
code2_code  = "C2X";

c = 299792458;
phase1_wavelength = c/(1575.42*10^6);
phase2_wavelength = c/(1227.60*10^6);
alfa = (1575.42*10^6)^2/(1227.60*10^6)^2;

n_epochs_to_display = nepochs;

% map container to map from GNSS system code to actual name
GNSSnames = containers.Map(["G","R","C","E"], ["GPS", "GLONASS", "Beidou", "Galileo"]);
nGNSSsystems = length(desiredGNSS_systems);

MP1 = []; MP2 = [];
SV_list = [];

for k = 1:nGNSSsystems
   %get current GNSS system
   sys = desiredGNSS_systems(k);
   
   % get current GNSS system index, return None if current sys is not in
   % data
   sysIndex = find([GNSSsystems{:}]==sys); 
   
   %fprintf('\nCurrently showing %s linear combinations\n\n', GNSSnames(sys))
   
   % iterate through epochs for the current GNSS system
   all_sat = NaN(n_epochs_to_display,32);
   fasebrudd = zeros(n_epochs_to_display,32);
   for epoch =1:n_epochs_to_display
      % get number of sat for current GNSS system that have observations
      % this epoch
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
         
         % get phaseobservations anf convert from units of cycles to meters
         phase1 = GNSS_obs{sysIndex}(SV, phase1_index, epoch)*phase1_wavelength;
         phase2 = GNSS_obs{sysIndex}(SV, phase2_index, epoch)*phase2_wavelength;
         
         % get code observations (m)
         code1 = GNSS_obs{sysIndex}(SV, code1_index, epoch);
         code2 = GNSS_obs{sysIndex}(SV, code2_index, epoch);
         
         %check that none of the phase obs are missing
         if ~any([phase1, phase2] == 0)
             
            % Fasebruddindikator
            IOD = (alfa/(alfa-1))*(phase1 - phase2); 
            % Multipath + Bias
            Mp1 = code1 - (1 + 2/(alfa-1))*phase1 + (2/(alfa-1))*phase2;
            Mp2 = code2 - (1 + 2/(alfa-1))*phase1 + (2/(alfa-1))*phase2;                     
           
         end
        
         all_sat(epoch,SV) = Mp1;
         fasebrudd(epoch, SV) = IOD; 
      end
      
   end
end


[m,n] = size(fasebrudd);
d = [];
index = [];
for i = 2:m
    for j = 1:n
        d(i,j) = fasebrudd(i,j);
        dis = d(i,j) - d(i-1, j);
        
        if abs(dis)  > 4/60
            index = [index;j,i];
                    
        end
    end
end  
index = sortrows(index,1);
[n, m] =size(index);

mean_table = [index(1,1),1,index(1,2),mean(all_sat(1:index(1,2),index(1,1)))];

for k = 1:n-1
    % Ser om det er mer enn ett fasebrudd, beregner kun gjennomsnitt mellom
    % fasebruudd
    if k > 1 && (index(k,1) == index(k+1,1)) && (index(k,1) ~= index(k-1,1))
        mean_table = [mean_table;index(k,1),1,index(k,2),mean(all_sat(1:index(k,2),index(k,1)),'omitnan')]; 
    
    end
    if (index(k,1) == index(k+1,1))  
        mean_table = [mean_table;index(k,1),index(k,2),index(k+1,2),mean(all_sat(index(k,2):index(k+1,2),index(k,1)),'omitnan')];
        %disp(index(k,1))
        %disp(index(k,2))
        %disp(index(k+1,2))
       % disp(mean(all_sat(index(k,2):index(k+1,2),index(k,1)),'omitnan'))
        %disp("--------------------------------------------------")
     % Om det kun er et fasebrud isf. Beregne det et gjennomsnitt for hele serien
    elseif (index(k,1) ~= index(k+1,1))
          mean_table = [mean_table;index(k,1),index(k,2),nepochs,mean(all_sat(index(k,2):nepochs,index(k,1)),'omitnan')]; 
    end
    
end
    
 
% Fjerner Bias fra MP (MP-Bias)    
[m,n] = size(mean_table);
real_MP = all_sat;
for i = 1:m
    real_MP(mean_table(i,2):mean_table(i,3),mean_table(i,1)) = abs(real_MP(mean_table(i,2):mean_table(i,3),mean_table(i,1))) - abs(mean_table(i,4));
end
    
% Fjerner alle kolonner med kun "NaN" verdier
real_MP = real_MP(:,~all(isnan(real_MP)));

% Plot
[m, n] = size(real_MP);

% Test
for i = 1:n
    plot(real_MP(:, i))
    hold on
end

%plot(real_MP(:,4))




