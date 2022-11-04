function type = frondtype(Nf,farm)
% Categorizes the frond into one of two categories: (1) SUBSURFACE = frond
% has not reachd the surface, criteria is that Nf at cultivation depth > 0
% (2) CANOPY = frond has reached the surface, criteria is that Nf at
% surface is > 0 
%
% Output:
%   type(subsurface,canopy,senescing)
%   1 = subsurface
%   2 = canopy
           
global param

              
%% SUBSURFACE                    
% if Nf at cultivation depth is > 0 = subsurface frond

   %DPD edit
   type(Nf(1) >0) = 1;
   %type(Nf(farm.z_arr == -farm.z_cult) >0) = 1;
   %type(Nf(farm.z_cult) > 0) = 1;

%% CANOPY
% if Nf at surface is > 0 = canopy frond; this will replace
% replaces subsurface category with canopy

   type(Nf(farm.z_arr > -farm.canopy) > 0) = 2;


end    
