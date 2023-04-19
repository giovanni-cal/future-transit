%% function that compute the 6 cost components and sums up them in " z "
function z = sum_local_costs(DX,x,r,R,v,v_p,tau_s,tau_p,rho_0,gamma,...
    c,mu,tot_Dem,ns,f)

        [z_L,z_V,z_M,z_A,z_W,z_T] = local_cost_fun(DX,x,r,R,v,v_p,...
            tau_s,tau_p,rho_0,gamma,c,mu,tot_Dem,ns,f);
        z = z_L + z_V + z_M + z_A + z_W + z_T;
        
end