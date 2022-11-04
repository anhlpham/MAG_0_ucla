function [UptakeN, UptakeFactor] = uptake(kelp,envt,farm,envt_counter)
% Determination of uptake rate for nitrate, ammonium, and urea.
%
%   Input: (ENVT,Q,Type,farm)
%           NO3 = nitrate concentration in seawater [umol NO3/m3]
%           NH4 = ammonium concentration in seawater [umol NH4/m3]
%           DON = dissolved organic nitrogen concentration [mmol/m3]
%           mag_u = seawater velocity [m/h]; magnitude velocity of x,y,z
%           Tw = wave period [h]
%           Type = 'subsurface/canopy' for conversion from surface area to
%           g(dry)
%           FARM_DIM = dimensions of farm needed for allometric
%           equations
%
%    Output:
%           Uptake, [mg N/g(dry)/h]; converted from [umol NO3+NH4+N/m2/h]
%           UptakeFactor.UptakeFactor_NO3/NH4/DON, dimensionless; 0-1 
%           UptakeFactor.Uptake_NOS/NH4/DON_mass, [umol NO3;NH4;N /g(dry)/h]
%
%    Note -> Uptake is only calculated for Canopy and Subsurface fronds. If
%    fronds are senscing, uptake not determined (NaN).


global param


%% ENVT INPUT

NO3 = envt.NO3(1:farm.nz,envt_counter);
NH4 = envt.NH4(1:farm.nz,envt_counter);
DON = envt.DON(1:farm.nz,envt_counter);
magu = envt.magu(1:farm.nz,envt_counter);
Tw = envt.Tw(1,envt_counter);


%% Q                                        
% Quota-limited uptake based on Droop "luxury uptake" formulation: 
% maximum uptake when nutrient quota Q is minimum and
% approaches zero as Q increases towards maximum; Possible that Q
% is greater than Qmax. Set any negative values to zero.

        UptakeFactor.vQ = (param.Qmax-kelp.Q)./(param.Qmax-param.Qmin);
        
        % Ensure that uptake doesn't take a negative value
        UptakeFactor.vQ(UptakeFactor.vQ < 0) = 0;
        UptakeFactor.vQ(UptakeFactor.vQ > 1) = 1;

        % Michaelis-Menten: Kinetically limited uptake only, for
        % exploratory, information purposes only. Not used to determine
        % Uptake rates
                
                 %vMichaelis_NO3 = NO3./(param.KsNO3+NO3); 
                 %vMichaelis_NH4 = NH4./(param.KsNH4+NH4);
                 %vMichaelis_DON = DON./(param.KsDON+DON);

%% v(Tw,Hs,[N])                 
% Kinetic+Mass-Transfer Limited Uptake: derivation from Stevens and
% Hurd 1997, Eq. 11
% The limitation of uptake due to kinetic and mass-transfer limitations was adapted from Stevens and Hurd (1997). 
% Essentially, it is a Michaelis-Menten parameterization of uptake that is depressed by slowed flow thickening 
% the diffusive boundary layer. Uptake is then replenished by wave action that acts by stripping away the diffusive boundary layer

        % CONSTANTS

            % visc = kinematic viscosity; [m2/h] 

            % Dm = molecular diffusion coefficient; [m2/h] *** Dm = T *
            % 3.65e-11 + 0.72e-10; determined for 18C (Li and Gregory
            % 1974);

            % n_length = number of iterations of wave-driven currents;
            % selected based on number of iterations attaining 95% of
            % maximum value for n_length = 1000

                % Viscosity in m2/hour
                visc     = 10^-6 * 3600;
                % Molecular diffusion coefficient in m2/hour
                Dm       = (18*3.65e-11 + 9.72e-10) * 3600;
                % Define number of terms for summation in equation for uptake
                n_length = 25; 

            % Diffusive Boundary Layer [m] Stevens and Hurd (1997) argued
            % that the high coefficient of u_star (0.33*u) was required to
            % match whole-frond drag estimates at low velocities

                DBL = 10 .* (visc ./ (0.33 .* abs(magu)));
                %disp('DBL'), DBL
        % There are two components of flow being considered. Calculate
        % contribution of each to uptake.

            % 1. Oscillatory Flow
                % Summation over n_length terms
                % DPD edit
                val = NaN(farm.nz,n_length);
                for n = 1:n_length
                    val(:,n) = (1-exp((-Dm.*n^2*pi^2.*Tw)./(2.*DBL.^2)))/(n^2*pi^2);
                end

                Oscillatory = ((4.*DBL)./Tw) .* sum(val,2);
                %disp('sum val'), sum(val,2)                
            % 2. Uni-directional Flow

                Flow = Dm ./ DBL;
		
        % Mass-Transfer Limitation is the sum of these two types of flows

                Beta   = Flow + Oscillatory; 
                %disp('Beta'), Beta
        % Kinetic+Mass-Transfer Limited Uptake, bring different components
        % to calcualte the Uptake Factor from Eq. 11 from Stevens and Hurd
        % 1997

                gammaNO3 = 1 +  (param.VmaxNO3 ./ (Beta.*param.KsNO3)) - (NO3 ./ param.KsNO3);
                gammaNH4 = 1 +  (param.VmaxNH4 ./ (Beta.*param.KsNH4)) - (NH4 ./ param.KsNH4);
                gammaDON = 1 +  (param.VmaxDON ./ (Beta.*param.KsDON)) - (DON.*0.2.*1e3 ./ param.KsDON); % DON*0.2 = [urea]
                
                % Below is what we call "Uptake Factor." It varies betwen 0
                % and 1 and includes kinetically limited uptake and
                % mass-transfer-limited uptake (oscillatory +
                % uni-directional flow)

                UptakeFactor.UptakeFactor_NO3 = NO3 ./ (param.KsNO3 .* ((NO3./param.KsNO3)  + 1/2 .* (gammaNO3+sqrt(gammaNO3.^2 + 4 .* (NO3 ./ param.KsNO3)))));
                UptakeFactor.UptakeFactor_NH4 = NH4 ./ (param.KsNH4 .* ((NH4./param.KsNH4)  + 1/2 .* (gammaNH4+sqrt(gammaNH4.^2 + 4 .* (NH4 ./ param.KsNH4)))));
                % DB: ### Check the conversion factors and units applied to DON ###
                UptakeFactor.UptakeFactor_DON = DON.*0.2.*1e3 ./ (param.KsDON .* ((DON.*0.2.*1e3./param.KsDON)  + 1/2 .* (gammaDON+sqrt(gammaDON.^2 + 4 .* (DON.*0.2.*1e3 ./ param.KsDON)))));
                
        % Uptake Rate CALCULATION, for each N source

                % Nutrient Uptake Rate = Max Uptake * UptakeFactor
                % [umol N/m2/h]

                    Uptake_NO3 = param.VmaxNO3 .* UptakeFactor.UptakeFactor_NO3; 
                    Uptake_NH4 = param.VmaxNH4 .* UptakeFactor.UptakeFactor_NH4; 
                    Uptake_DON = param.VmaxDON .* UptakeFactor.UptakeFactor_DON; 

                % Convert from surface area to g(dry). Based on
                % allometric conversions from param that
                % are dependent on kelp type (subsurface, canopy,
                % watercolumn)
                % [umol N/g(dry)/h]

                % and multiply by Q limitation (vQ)
                    
                    UptakeFactor.Uptake_NO3_mass = Uptake_NO3 .* kelp.sa2b.*2 ./ param.dry_wet .* UptakeFactor.vQ;
                    UptakeFactor.Uptake_NH4_mass = Uptake_NH4 .* kelp.sa2b.*2 ./ param.dry_wet .* UptakeFactor.vQ;
                    UptakeFactor.Uptake_DON_mass = Uptake_DON .* kelp.sa2b.*2 ./ param.dry_wet .* UptakeFactor.vQ;
               
                
                % Convert from umol -> mg N
                % [mg N/g(dry)/h]

                    Uptake_NO3_massN = UptakeFactor.Uptake_NO3_mass .* param.MW_N ./ 1e3; 
                    Uptake_NH4_massN = UptakeFactor.Uptake_NH4_mass .* param.MW_N ./ 1e3;
                    Uptake_DON_massN = UptakeFactor.Uptake_DON_mass .* param.MW_N ./ 1e3;


%% TOTAL Uptake = Uptake NO3 + Uptake NH4 + Uptake DON
% [mg N/g(dry)/h]

    UptakeN = Uptake_NO3_massN + Uptake_NH4_massN + Uptake_DON_massN;
    UptakeN(kelp.Nf == 0) = 0;
   
end
