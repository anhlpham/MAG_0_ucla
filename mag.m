function [kelp, DON, PON] = mag(kelp,envt,farm,time,envt_counter)
% Calculate uptake and growth to derive delta_Ns and delta_Nt
% function dependencies: uptake, growth
%
% Input: kelp_fr, envt, farm, time, envt_step, growth_step
% Output: kelp_fr at t+1


global param
                
%-------------------------------------------------------------------------
% KELP - Known Variables
% create temporary variables

    Ns = kelp.Ns;
    Nf = kelp.Nf;
                
%-------------------------------------------------------------------------
% UPTAKE
% [mg N/g(dry)/h]: NO3+NH4+Urea

    Uptake = uptake(kelp,envt,farm,envt_counter);   
    
%-------------------------------------------------------------------------
% GROWTH
% Growth limited by internal nutrients, temperature, PAR, and height (stops
% growing once reaches max height); [h-1]

    Growth = growth(kelp,envt,farm,envt_counter);       
           
%-------------------------------------------------------------------------
% MORTALTIY                
% d_wave = frond loss due to waves; dependent on Hs, significant
% wave height [m]; Rodrigues et al. 2018 demonstrates linear relationship
% between Hs and frond loss rate [h-1] (continuous)

    M_wave  = param.d_wave_m .* envt.Hs(1,envt_counter);

% d_blade = blade erosion [h-1] (continuous); Multiplied by the
% frBlade -> fraction of total as blade
    
    M_blade = param.d_blade .* kelp.frBlade;  
    M_tot = M_wave + M_blade;
    
% Ns at t+1
% Because Ns is redistributed based on distribution of Nf -> we calculated
% Nf(t+1) first and now follow-up with Ns

% Ns(t+1) = Ns(t) + Uptake - Growth - Mortality 
% For uptake, only biomass in blades (frBlade) contributes 
    %disp('frBlade'), kelp.frBlade
    dNs1 = Uptake .* kelp.B .* kelp.frBlade .* time.dt_Gr ... % uptake contribution from blades only
         - Growth .* Ns .* time.dt_Gr ... % stored nitrogen lost due to growth
         - param.d_dissolved .* Ns .* time.dt_Gr ... % exudation term
         - M_tot .* Ns .* time.dt_Gr; % wave-based mortality 
    Ns_new = Ns + dNs1;  % add dNs to Ns

% Nf at t+1
            
% Nf(t+1) = Nf(t) + Growth - Mortality

     % change in Nf
     dNf1 = Growth .* Ns .* time.dt_Gr ...
          - M_tot .* Nf .* time.dt_Gr;
     Nf_new = Nf + dNf1; % add dNf to Nf

     
%-------------------------------------------------------------------------
% VERTICAL DISTRIBUTION OF STATE VARIABLES
% Nf distributed based on height-biomass relationships derived from
% mag1_frond and saved in a table, b_per_m.mat. This table has been loaded
% to param.b_per_m

    temp_Nf = find_nan(Nf_new);
    Nf_new = trapz(farm.z_arr,temp_Nf) .* kelp.b_per_m;
    
    % Ns redistributed same as Nf to maintain Q along frond (translocation)
    
    temp_Ns = find_nan(Ns_new);
    Ns_new = trapz(farm.z_arr,temp_Ns) .* kelp.b_per_m;



%-------------------------------------------------------------------------
% DON, PON

% % DON, the amount of kelp contribution to DON; [mg N] -> [mmol N/m3]
% %   There are three sources of DON contribution from the Ns pool
% %   1. dissolved loss
% %   2. blade erosion
% %   3. wave-based mortality
% 
%     DON = (nansum(param.d_dissolved .* Ns) .* time.dt_Gr) ...
%         +  nansum(M_tot .* Ns) .* time.dt_Gr ...
%         ./ param.MW_N;
%   
%    
% % PON, the amount of kelp contribution to PON; [mg N/m3]
% %   There are three sources of PON contribution from the Nf pool
% %   1. blade erosion
% %   2. wave-based mortality
% 
%     PON = nansum(M_tot .* Nf) .* time.dt_Gr;
                          
%-------------------------------------------------------------------------
% UPDATE STATE VARIABLES
% Note, these replace existing matrices (don't append)

kelp.Nf = Nf_new;
kelp.Ns = Ns_new;

%disp('plant height'), kelp.height

end
