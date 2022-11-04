function kelp = frondsenescence(kelp,time,sim_hour)
% Identify those fronds that are senescing. This is defined as 2 weeks past
% maximum age. Nf and Ns is lost at a rate of 10% for 10 days for the
% fraction of total fronds that are senescing

global param

    % assess frond status relative to age
    kelp.fronds.status(sim_hour < kelp.fronds.end_age) = 1; % alive
    kelp.fronds.status(sim_hour > kelp.fronds.end_age) = 2; % senesce
    kelp.fronds.status(sim_hour > kelp.fronds.end_age+24*10) = 3; % dead
    
    % fraction of senescing fronds
    s = sum(kelp.fronds.status==2);
    a = sum(kelp.fronds.status==1);
    f_s = s/(s+a);
    disp('Time'), sim_hour / 24
    disp('s'), s
    disp('a'), a
    disp('f_s'), f_s
    
    % remove Nf and Ns at a rate of d_frond for fraction that are senescing
    Nf_senesce = param.d_frond .* kelp.Nf .* f_s .* time.dt_Gr;
    Ns_senesce = param.d_frond .* kelp.Ns .* f_s .* time.dt_Gr;
    
    % update state variables
    kelp.Nf = kelp.Nf - Nf_senesce;
    kelp.Ns = kelp.Ns - Ns_senesce;
    
end
