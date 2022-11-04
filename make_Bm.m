function [Bm_out] = make_Bm(h_p,farm) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create analytical form of biomass per meter (Bm)
% that are normalized so they integrate to 1
% Bm(z) depends on kelp height and a canopy threshold
% h_p : height of the plant
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global param
% Declare an empty b_per_m array
Bm = NaN(farm.nz,1);
alpha = 2; %scale of vertical transition from subsurface to canopy
% Get depth of top of the plant in water column
z_p = -farm.z_cult + h_p;
%disp('z_p'), z_p
% Check if z_p exceeds surface
if z_p>0
   z_p = 0;
end

%Create Bm
%Subsurface Case
if z_p < param.z_canopy
    %sgn_p = ((sign(z_p - farm.z_arr+farm.dz) + 1) / 2.);
    %Bm(:) = sgn_p * farm.z_arr^2   ; 
       for k=1:length(Bm)
       %this first if is not necessary in matlab code 
       % b/c depths do not go below farm itself
       if farm.z_arr(k)<-farm.z_cult
	  Bm(k) = 0;
       %In farm, below plant height
       elseif farm.z_arr(k)<=z_p
	  Bm(k) = farm.z_arr(k)^2;
       %In farm, above plant height
       else
	  Bm(k) = 0;
       end
   end
end

%Canopy Case
if z_p>=param.z_canopy
   %Bm(:) = (1./abs(farm.z_cult)) + exp((farm.z_arr - param.z_canopy) / vscale);
   for k=1:length(Bm)
       %this first if is not necessary in matlab code 
       % b/c depths do not go below farm itself
       if farm.z_arr(k)<-farm.z_cult
	  Bm(k) = 0;
       else
	   Bm(k) = (1 / farm.z_cult) + exp((farm.z_arr(k) - param.z_canopy)/alpha);
       end
       %Not in canopy
       %elseif farm.z_arr(k)<=param.z_canopy
	  %Bm(k) = param.B0;
       %In canopy
       %else
%	  Bm(k) = exp(farm.z_arr(k));
 %      end
   end
end

%Normalize Bm so that it integrates to 1
int_B = trapz(farm.z_arr,Bm);
%disp('int_Bm'), int_B
Bm_out = Bm / int_B;


%Make constant for testing
%Bm_temp = NaN(farm.nz,1);
%Bm_temp(:) = 0.01;
%int_B = trapz(farm.z_arr,Bm_temp);
%Bm_out = Bm_temp / int_B;


