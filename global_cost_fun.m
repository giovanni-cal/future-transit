function [F_V,F_M,F_W,F_T,F_e_T] = global_cost_fun(x,r,phi_B,H_B)

% Commercial speed on boundary ring
    v_comm_B = 1 / (1 / v(1) + (tau_s(1) / (phi_B * x))); %or Theta_B ???
% DEMAND(x)
    Prob_x = @(x) Dem/tot_Dem;
    Prob_rR = integral(Prob_x,r,R);
    %% AGENCY COSTS
% Vehicle-distance traveled per hour: F_V
    F_V = (4 * pi * x) / H_B;
% Fleet size F_M
    F_M = (4 * pi * x) / (H_B * v_comm_B);

    %% USER COSTS
% Average waiting time: F_W
    F_W = Dem * 2/PI * H_B/2 * Prob_rR^2;
% In-vehicle travel time: F_T
    F_T = Dem * 2*r/pi * 1/v_comm_B * Prob_rR^2;
% number of transfers F_e_T
    F_e_T = Dem * (1 + (2/pi * Prob_rR^2));
end