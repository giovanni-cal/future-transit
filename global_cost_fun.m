function [F_V,F_M,F_W,F_T,F_transf] = global_cost_fun(r,phi_B,H_B,v,tau_s,c_V,c_M,mu_W,mu_A,mu_T,Prob_rR,tot_Dem)

% Commercial speed on boundary ring
    v_comm_B = 1 / (1 / v(1) + (tau_s(1) / (phi_B * r))); %or Theta_B ???
    
%capacity constraint in the "master" script

    %% AGENCY COSTS
% Vehicle-distance traveled per hour: F_V
    F_V = c_V(1) * (4*pi*r)/H_B;
% Fleet size F_M
    F_M = c_M(1) * (4*pi*r)/(H_B*v_comm_B);

    %% USER COSTS
% Average waiting time: F_W
    F_W = mu_W * tot_Dem * 2/pi * H_B/2 * Prob_rR^2;
% In-vehicle travel time: F_T
    F_T = mu_T * tot_Dem * 2*r/pi * 1/v_comm_B * Prob_rR^2;
% number of transfers F_transf
    penalty = 0.05; %3 minutes walk per transfer
    F_transf = mu_A * penalty * tot_Dem * (1 + (2/pi * Prob_rR^2));

end