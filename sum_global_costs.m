function F = sum_global_costs(r,phi_B,H_B,Th_r_min,v,tau_s,c,mu,Prob_rR,tot_Dem)

    [F_V,F_M,F_W,F_T,F_transf] = global_cost_fun(r,phi_B,H_B,Th_r_min,v,tau_s,...
        c,mu,Prob_rR,tot_Dem);

    F = F_V + F_M + F_W + F_T + F_transf;
    
end