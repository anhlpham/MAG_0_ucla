function kelp = kelpchar(kelp,farm)
% Calculate biological characteristics from Nf, Ns (known)
%
% OUTPUT:  
%   Q (## UNITS ###) mg N / g-dry of biomass
%   Biomass (## UNITS ###) g-dry of biomass
%   depth resolved surface_area-to-biomass (## UNITS ###) sa2b [m2/g(wet)]
%   height: total height (m)
%   frBlade, fractional biomass that is blade = blade:frond (unitless)
%
% NOTES:
% Q is integrated across depth, Rassweiler et al. (2018) found that N content does
% not vary with depth. This can be interpreted as translocation occuring on
% the scale of hours (Parker 1965, 1966, Schmitz and Lobban 1976). 
% Side note: Ns is redistributed after uptake as a function of fractional Nf.
% This is how the model "translocates" along a gradient of high-to-low
% uptake. Mathematically, this keeps Q constant with depth. -> There may be
% more recent work coming out of SBC LTER that indicates N varies along
% the frond, particularly in the canopy in a predictable age-like matter
% (e.g., what tissues are doing most of the photosynthesis) (T. Bell pers.
% comm.)
%
% Biomass calculation from Hadley et al. (2015), Table 4 [g(dry)/dz]
%                
% Surface area to biomass
%
% Blade-to-stipe ratio derived from Nyman et al. 1993 Table 2
  

global param
% Calculate DERIVED variables on a per frond basis

    % KNOWN STATE VARIABLES
    % Ns, Nf 
    
    % DERIVED VARIABLES
        
        % First set NaN values to zeros in the temporary arrays
        temp_Ns = find_nan(kelp.Ns);
        temp_Nf = find_nan(kelp.Nf);
        
        % Calculates the nutrient quota Q (excess stored nutrients, unitless)
        kelp.Q = param.Qmin .* (1 + trapz(farm.z_arr,temp_Ns) ./ trapz(farm.z_arr,temp_Nf));
        kelp.B = kelp.Nf ./ param.Qmin; % grams-dry m-3
        
        % Defines surface-area to biomass conversion (depth resolved)
        % sa2b [m2/g(wet)]
        % Note: canopy forming has more surface area in canopy
        % Two cases:
        % (1) Canopy is present: uses two conversion factors:
	%     - surface area of canopy
        %     - surface area of plant below the canopy
        % (2) Canopy not present, all plant submerged
        %     - surface area of submerged plant
        kelp.sa2b = NaN(farm.nz,1);
        if any(kelp.Nf(farm.z_arr > -farm.canopy) > 0)
            kelp.sa2b(farm.z_arr >= -farm.canopy) = param.Biomass_surfacearea_canopy;
            kelp.sa2b(farm.z_arr < -farm.canopy) = param.Biomass_surfacearea_watercolumn;
        else
            kelp.sa2b(1:farm.nz) = param.Biomass_surfacearea_subsurface;
        end
        
        % First set NaN values to zeros in the temporary array temp_B
        temp_B = find_nan(kelp.B);
        kelp.height = (param.Hmax .* trapz(farm.z_arr,temp_B)./1e3 )./ (param.Kh + trapz(farm.z_arr,temp_B)./1e3);
        
        % Blade to Stipe for blade-specific parameters
           
            % generate a fractional height
            fh = diff(farm.z_arr)';
            fh(end+1) = fh(end);
            fh = cumsum(fh);
            fh = fh .* (kelp.B > 0);
            %fh = fh .* ~isnan(kelp.B);
            %disp('fh pre'), fh
            fh = fh ./ kelp.height; 
            fh(fh==0) = NaN; 
            fh(fh>1) = 1;
	    BtoS = param.Blade_stipe(1) - param.Blade_stipe(2) .* fh + param.Blade_stipe(3) .* fh .^ 2;
	    
            kelp.frBlade = BtoS ./ (BtoS + 1);
            clear fh BtoS
            
        % biomass per m (for growth)
        % (shape function (## UNITS ###)
        
        kelp.b_per_m = make_Bm(kelp.height,farm);
    
end
