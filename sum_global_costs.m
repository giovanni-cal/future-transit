function F = sum_global_costs(r,phi_B,H_B,v,tau_s,c_V,c_M,mu_W,mu_A,mu_T,Prob_rR,tot_Dem)

    [F_V,F_M,F_W,F_T,F_transf] = global_cost_fun(r,phi_B,H_B,v,tau_s,...
        c_V,c_M,mu_W,mu_A,mu_T,Prob_rR,tot_Dem);

    F = F_V + F_M + F_W + F_T + F_transf;
    
end