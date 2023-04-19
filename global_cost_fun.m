function [F_V,F_M,F_W,F_T,F_transf] = ...
    global_cost_fun(r,phi_B,H_B,Th_r_min,v,tau_s,c,mu,Prob_rR,tot_Dem)

% Commercial speed on boundary ring
    v_comm_B = 1 / (1 / v(1) + (tau_s(1) / (phi_B * r))); %or Theta_B ???
    
%capacity constraint in the "master" script

    %% AGENCY COSTS
% Vehicle-distance traveled per hour: F_V
    F_V = c(2,1) * (4*pi*r)/H_B;
% Fleet size F_M
    F_M = c(3,1) * (4*pi*r)/(H_B*v_comm_B);

    %% USER COSTS
% Average waiting time: F_W
    F_W = mu(2,1) * tot_Dem * (2-Th_r_min)/pi * H_B/2 * Prob_rR^2;
% In-vehicle travel time: F_T
    F_T = mu(1,1) * tot_Dem * (2-Th_r_min)*r/pi * 1/v_comm_B * Prob_rR^2;
% number of transfers F_transf
    penalty = 0.033; %2 minutes walk per transfer
    F_transf = mu(3,1) * penalty * tot_Dem * (1 + ((2-Th_r_min)/pi * Prob_rR^2));

end