function envt = envt_sb(farm,time,dir_ROMS,dir_WAVE)
% Environmental Input Data
% INPUT: farm, time, directories for ROMS and WAVE data
% OUTPUT: envt.()
%   NO3; daily ROMS Nitrate in seawater; [umol/m3]
%   NH4; daily ROMS Ammonium in seawater; [umol/m3]
%   DON; daily ROMS DON, dissolved organic nitrogen; [mmol N/m3]
%   T; daily ROMS Temperature; [Celcius]
%   magu; magnitude velocity from daily ROMS Uo,Vo,Wo Seawater velocity, [m/h]
%   PAR; daily ROMS PAR, photosynthetically active radiation; incoming PAR; [W/m2]
%   chla; adily ROMS chl-a: sum of DIAZ+DIAT+SP (small phytoplankton); [mg-chla/m3]
%   Tw; daily NDCP Wave period; [h]
%   Hs; daily NDCP Significant wave height; [m]
%
% The ROMS simulations have been pre-processed and saved as mat files for
% use by MAG.
% Data is 1-m bins from surface to bottom depth.
% Interpolate to the farm z-array.
%
% The NDCB data has been downloaded and saved as mat files for use by MAG.

%% File Directory
           
% SBC farm site; offshore Mohawk; 60-m water depth

    ROMS_start = datenum([1994 01 01 0 0 0]); % from ROMS folks
    
% Extract days

    filename = strcat(dir_ROMS,'NO3.mat'); NO3 = load(filename);
    ROMS_time = datevec(ROMS_start + NO3.ocean_time ./ (60*60*24)); clear NO3
    ROMS_time = datenum(ROMS_time(:,1:3));
    
    start = min(time.timevec_ROMS);
    stop  = max(time.timevec_ROMS);
    [~,idx_start]=min(abs(ROMS_time-start));
    [~,idx_stop] =min(abs(ROMS_time-stop));
    
    % the date indices to extract from full ROMS time series
    ROMS_extract  = idx_start:idx_stop;
    clear val start stop filename ROMS_time idx_start idx_stop
    
    % depth domain of ROMS
    ROMS_z = 0:1:59;
    
% Load ROMS data and create variable to store boundary condition.
%% Nitrate 
            
    filename = strcat(dir_ROMS,'NO3.mat');
    NO3 = load(filename);

    NO3 = NO3.NO3(ROMS_extract,:);
    for rr = 1:length(ROMS_extract)
        envt.NO3(rr,:) = interp1(ROMS_z,NO3(rr,:),abs(farm.z_arr));
    end
    envt.NO3 = envt.NO3' .* 1e3; % umol/m3
    envt.NO3(envt.NO3 <= 0.01e3) = 0.01e3; % replace negatives
    
        clear NO3 filename rr

        
%% Ammonium
        
    filename = strcat(dir_ROMS,'NH4.mat');
    NH4 = load(filename);

    NH4 = NH4.NH4(ROMS_extract,:);
    for rr = 1:length(ROMS_extract)
        envt.NH4(rr,:) = interp1(ROMS_z,NH4(rr,:),abs(farm.z_arr));
    end
    envt.NH4 = envt.NH4' .* 1e3; % umol/m3
    
        clear NH4 filename rr

        
%% DON
        
    filename = strcat(dir_ROMS,'DON.mat');
    DON = load(filename);
    
    DON = DON.DON(ROMS_extract,:);
    for rr = 1:length(ROMS_extract)
        envt.DON(rr,:) = interp1(ROMS_z,DON(rr,:),abs(farm.z_arr));
    end
    envt.DON = envt.DON';
    
        clear DON filename rr


%% Temperature
    
    filename = strcat(dir_ROMS,'temp.mat');
    temp = load(filename);
    
    T = temp.temp(ROMS_extract,:);
    for rr = 1:length(ROMS_extract)
        envt.T(rr,:) = interp1(ROMS_z,T(rr,:),abs(farm.z_arr));
    end
    envt.T = envt.T';
    
        clear temp filename T

            
%% Seawater Velocity, u,v,w
% Seawater velocity in x
            
    filename = strcat(dir_ROMS,'u.mat');
    u = load(filename);
    
    u= u.u(ROMS_extract,:);
    u = u .* 60 .* 60; % [m/h]
    
        clear filename

% Seawater velocity in y
            
    filename = strcat(dir_ROMS,'v.mat');
    v = load(filename);
    
    v = v.v(ROMS_extract,:);
    v = v .* 60 .* 60; % [m/h]
    
        clear filename
            
% Seawater velocity in z
% This is informed by ROMS and based on conversation with Kristen on
% 20191015 going to set vertical velocities to zero

    filename = strcat(dir_ROMS,'w.mat');
    w = load(filename);
    
    w = w.w(ROMS_extract,:);
    w = w .* 60 .* 60; % [m/h]
    w = w .* 0;
    
        clear filename

% Seawater magnitude velocity        
    magu = sqrt(u.^2 + v.^2 + w.^2);
    for rr = 1:length(ROMS_extract)
        envt.magu(rr,:) = interp1(ROMS_z,magu(rr,:),abs(farm.z_arr));
    end
    envt.magu = envt.magu';
    
clear u v w rr magu
    
%% PAR

    % Load ROMS data of PAR and extract surface value only
    % Also load PARincoming which isn't all that different from PAR at
    % the surface. PARincoming is 0.45 * penetration of solar heat
    % PAR at surface is modified as = PARinc * (1 - exp(Kpar))
    % [W/m2]
        
    filename = strcat(dir_ROMS,'PAR.mat');
    PAR = load(filename);
    
    envt.PAR = PAR.PAR(ROMS_extract,1)';
    
    
        clear PAR filename
           
            
%% CHL-a
% sum of three phytoplankton components

    filename1 = strcat(dir_ROMS,'DIATCHL.mat');
    filename2 = strcat(dir_ROMS,'DIAZCHL.mat');
    filename3 = strcat(dir_ROMS,'SPCHL.mat');
    DIATCHL = load(filename1);
    DIAZCHL = load(filename2);
    SPCHL = load(filename3);

    chla = ...
          DIATCHL.DIATCHL(ROMS_extract,:)...
        + DIAZCHL.DIAZCHL(ROMS_extract,:)...
        + SPCHL.SPCHL(ROMS_extract,:);
    chla(chla < 0) = 0; % replace negatives
    for rr = 1:length(ROMS_extract)
        envt.chla(rr,:) = interp1(ROMS_z,chla(rr,:),abs(farm.z_arr));
    end
    envt.chla = envt.chla';
    
        clear DIATCHL DIAZCHL SPCHL filename1 filename2 filename3 chla rr
            
            
%% Wave period, Significant wave height
NDCPfilename = strcat(dir_WAVE,'NDCP46053.mat');

    % NDCP; mat file
    % Column 1 = matlab time vec
    % Column 2 = Hs
    % Column 3 = Tw
    wave = load(NDCPfilename);
     
    envt.Tw = interp1(wave.NDCP46053(:,1),wave.NDCP46053(:,3),time.timevec_ROMS);
    envt.Tw = envt.Tw ./ (60*60); % [h]
    envt.Hs = interp1(wave.NDCP46053(:,1),wave.NDCP46053(:,2),time.timevec_ROMS);

       clear wave
         
clear ROMS_dir ROMS_start ROMS_extract 
end