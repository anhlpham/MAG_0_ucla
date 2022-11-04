function [Growth, gQ, gT, gE, gH] = growth(kelp,envt,farm,envt_counter)
% Growth, nitrogen movement from Ns to Nf 
% Growth = umax * gQ * gT * gE * gH; [per hour]
%
% Input: (Q,B,envt,farm,envt_counter)
%
% Output: 
%   Growth, [h-1]
%   gQ, quota-limited growth 
%       from Wheeler and North 1980 Fig. 2   
%   gT, temperature-limited growth
%       piecewise approach taken from Broch and Slagstad 2012 (for sugar
%       kelp) and optimized for Macrocystis pyrifera
%       calls on gT_vX
%   gE, light-limited growth
%       from Dean and Jacobsen 1984
%       calls on gE_vX
%   gH, carrying-capacity limited growth, from WASP manual Equation 33



global param
%% gQ -> ranges from zero to 1

    gQ = (kelp.Q - param.Qmin) ./ (param.Qmax - param.Qmin); 
    gQ(gQ > 1) = 1;
    gQ(gQ < 0) = 0;
    
    
%% gT -> ranges froms zero to 1

    % temp data
    temp = envt.T(1:farm.nz,envt_counter);
    
    gT = NaN(size(temp)); % preallocate space
        
        gT(temp < param.Tmin) = 1/param.Tmin * temp(temp < param.Tmin);
        gT(temp >= param.Tmin & temp < param.Tmax) = 1;
        
            % Solve systems of equations where Tmax = 1; Tlim = 0;
            % Linear decrease from Tmax to Tlim with intercept b and slope m
            
            b = 1 / (-param.Tmax / param.Tlim + 1);
            m = 1 / (param.Tmax - param.Tlim);
            %disp('b'), b
	    %disp('m'), m

        gT(temp >= param.Tmax & temp <= param.Tlim) = m .* temp(temp >= param.Tmax & temp <= param.Tlim) + b;
        gT(temp > param.Tlim) = 0;

    clear temp
    
    
%% gE -> ranges from zero to 1
% light varies across the farm, so extract correct light field with kelploc
% Bertalanffy Growth Equation (Dean and Jacobsen 1984)
% Input: PAR [W/m2]

    % k sets curvature of relationship. qualitatively fit to match Dean
    % and Jacobsen 1984; 50% growth at PAR ~2.5 and near 100% growth at
    % PAR ~7+
        
        gE = 1-exp(-param.kPAR*(envt.PARz-param.PARc));
        %daniel hardocde for testing
	%gE = 1 - exp(param.kPAR*(86-param.PARc))
        % If values < 0 replace with zero. We are explicitely modeling
        % mortality and so growth shouldn't be negative.
        gE(gE < 0) = 0;
	%disp('PAR'), envt.PARz
       % disp('gE'), gE      
        
%% gH -> ranges from zero to 1           
% as frond approaches carrying capacity; growth is limited (approaches
% zero) -> space limitation effect

    %gH = 1-(nansum(kelp.B)./param.kcap).^2; % in units of g-dry

    % Hyperbolic tangent function. Sets gH = 0.5 at height = 95% of max height
    gH = 0.5 + 0.5 .* tanh(-(kelp.height - 0.95 * param.Hmax));
    %disp('tanh height'), kelp.height
    %disp('matlab tanh'), tanh(-(kelp.height - (param.Hmax-0.05*param.Hmax)))
    %disp('gH'), gH
    %disp('gQ'), gQ
    %disp('gT'), gT
    %disp('gE'), gE
    %figure
    %plot(gE,farm.z_arr)
    %title('gE small dz')
    %plot(gT,farm.z_arr)
    %title('gT small dz')
    %plot(gH,farm.z_arr)
    %title('gH small dz')


%% Growth
% per hour

    %daniel testing absence of gE
    %Growth = param.umax .* gQ .* gT .* 1. * gH;
    Growth = param.umax .* gQ .* gT .* gE .* gH;
    Growth(isnan(kelp.Nf)) = NaN; % only retain values where kelp is present

end
