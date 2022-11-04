function time = simtime(duration)
% simulation time and intervals; unit is hours
% Input: duration as a vector [start; stop]
% Output: structure containing
%   start = datevec of simulation start
%   duration of simulation (in hours)
%   dt_Gr = time step for growth (in hours)
%   dt_ROMS = time step of environmental input
%   timevecus in matlab datenum format

    time.start      = duration(1,:);
    time.duration   = (datenum(duration(2,:))-datenum(duration(1,:)))*24; % days * hours/day
    time.dt_Gr      = 24; % solve growth every X hours
    %time.dt_Gr = 1;
    time.dt_ROMS    = 24; % based on extracted ROMS files
     
    time.timevec_Gr = datenum(time.start):time.dt_Gr/24:datenum(time.start)+time.duration/24-1/24;
    time.timevec_ROMS = datenum(time.start):time.dt_ROMS/24:datenum(time.start)+time.duration/24-1/24;
    
end  
