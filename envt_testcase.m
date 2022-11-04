function envt = envt_testcase(farm,time)
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
% The NDCB data has been downloaded and saved as mat files for use by MAG.

%% File Directory
             
    envt.NO3 = 1.557 .* 1e3; % umol/m3
    envt.NO3 = repmat(envt.NO3,farm.nz,length(time.timevec_Gr));
        
%% Ammonium
        
    envt.NH4 = 87.3; % umol/m3
    envt.NH4 = repmat(envt.NH4,farm.nz,length(time.timevec_Gr));
    
        
%% DON
        
    envt.DON = 3.07; % mmol/m3
    envt.DON = repmat(envt.DON,farm.nz,length(time.timevec_Gr));
   

%% Temperature
    
    envt.T = 14.5; % umol/m3
    envt.T = repmat(envt.T,farm.nz,length(time.timevec_Gr));
   
            
%% Seawater Velocity, u,v,w
% Seawater magnitude velocity    

    envt.magu = 516; % m/h
    envt.magu = repmat(envt.magu,farm.nz,length(time.timevec_Gr));
   

%% PAR

    envt.PAR = 86; % [W/m2]
    envt.PAR = repmat(envt.PAR,farm.nz,length(time.timevec_Gr));
           
            
%% CHL-a
% sum of three phytoplankton components

    envt.chla = 0.4; % 
    envt.chla = repmat(envt.chla,farm.nz,length(time.timevec_Gr));
          
            
%% Wave period, Significant wave height

    envt.Tw = 0.0021; % [h]
    envt.Tw = repmat(envt.Tw,farm.nz,length(time.timevec_Gr));
    
    envt.Hs = 0.86; % [m]
    envt.Hs = repmat(envt.Hs,farm.nz,length(time.timevec_Gr));
       
    
end
