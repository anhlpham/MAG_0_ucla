function farm = farmdesign
% structure with farm properties
%
% Output: farm.()
%   z_cult; depth of cultivation, [m]
%   x,y,z; dimension of farm [m]
%   dx,dy,dz; bin size [m] 
%   seeding; initial seeding biomass

%% farm dimensions and grid

    farm.z_cult = 20; %[m below surface]
    farm.z      = farm.z_cult; % [m]
    farm.z_arr = linspace(-farm.z,0,100);
    farm.dz = farm.z_arr(2) - farm.z_arr(1);
    farm.nz = length(farm.z_arr);
    
    % initial B/Q conditions
    farm.seedingB = 3*1e3; % seeding biomass [g-dry m-2] ## CHECK UNITS ### (it said per meter)
    farm.seedingQ = 15; % seeding Q
    
    % 'canopy' starts at what depth below the surface
    farm.canopy = 1; % what depth is canopy defined at... (m)
    
end
