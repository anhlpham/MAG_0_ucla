function PARz = canopyshading(kelp,envt,farm,envt_counter)
% Bio-optical model, calculates the attenuation of light due to water,
% chl-a (from ROMS) and nitrogen-specific shading. Incoming PAR from ROMS.
% 
% Input: Nf (per m3; not per frond) * already smoothed at canopy
%        envt, farm, step (ENVT data)
% 
% Output: PAR as a function of depth across the entire farm regardless of
% whether there is kelp present in a given area or not.
%
% DB 10/28/22: Check attenuation of light profile
%              seems light attenuate *faster* below the canopy -- should it be the other
%              way around?

global param

% Canopy Shading; Nf

% canopy shading is in fixed nitrogen units
Nf = kelp.Nf;

% redistribute amount of Nf at surface 
canopyHeight = kelp.height - farm.z_cult; 
canopyHeight(canopyHeight<1) = 1;
% disp('canopy height'), canopyHeight

% below I am defining "canopy" as depths less than 1 meter; ...
for k=1:length(farm.z_arr)
    if farm.z_arr(k) > -farm.canopy
        Nf(k) = Nf(k)  / canopyHeight; 
    end

    % Replacement of NaN with zero for mathematical reasons   
    Nf(isnan(Nf)) = 0; 

    % Attenuation of PAR with depth
    % PAR, incoming
    PAR0 = envt.PAR(1,envt_counter);

    % preallocate space
    PARz = NaN(farm.nz,1);

   % Calculate attenuation coefficents and resulting PAR from surface to
   % cultivation depth

    for indz = length(farm.z_arr):-1:1
        z = farm.z_arr(indz);
        if z == 0 
            PARz(indz) = PAR0; % no attenuation at surface
        else
            
        % attenuate with sum of three contributions
        K = param.PAR_Ksw + ...
            param.PAR_Kchla * envt.chla(indz,envt_counter) + ... 
            param.PAR_KNf * Nf(indz+1);

        PARz(indz) = PARz(indz+1) .* (exp(-K * farm.dz)); % output variable
               
        end
    end
end
