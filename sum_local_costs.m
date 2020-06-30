function z = sum_local_costs(DX,x,r,R,v,v_p,tau_s,tau_p,d,rho_0,gamma,...
    c_L,c_V,c_M,mu_W,mu_A,mu_T,tot_Dem,f)

    [z_L,z_V,z_M,z_A,z_W,z_T] = local_cost_fun(DX,x,r,R,v,v_p,...
        tau_s,tau_p,d,rho_0,gamma,c_L,c_V,c_M,mu_W,mu_A,mu_T,tot_Dem,f);

    z = z_L + z_V + z_M + z_A + z_W + z_T;
    
end