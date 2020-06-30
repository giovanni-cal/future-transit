function z = sum_local_costs(DX,x,r,R,v,v_p,tau_s,tau_p,d,rho_0,gamma,...
    c_L,c_V,c_M,mu_W,mu_A,mu_T,tot_Dem,f)

   if f== 2 || f==3 or ....
        [z_L,z_V,z_M,z_A,z_W,z_T] = local_cost_fun(DX,x,r,R,v,v_p,...
        tau_s,tau_p,d,rho_0,gamma,c_L,c_V,c_M,mu_W,mu_A,mu_T,tot_Dem,f);
   elseif f==6
        [z_L_fixed,z_V_fixed,z_M_fixed,z_A_fixed,z_W_fixed,z_T_fixed] = 
             local_cost_fun(DX,x,r,R,v,v_p,...
    tau_s,tau_p,d,rho_0,gamma,c_L,c_V,c_M,mu_W,mu_A,mu_T,tot_Dem,3)
    z_fixed = z_L_fixed + z_V_fixed + z_M_fixed + z_A_fixed + z_W_fixed + z_T_fixed;
    
        [z_L_demresp,z_V_demresp,z_M_demresp,z_A_demresp,z_W_demresp,z_T_demresp] = 
             local_cost_fun(DX,x,r,R,v,v_p,...
    tau_s,tau_p,d,rho_0,gamma,c_L,c_V,c_M,mu_W,mu_A,mu_T,tot_Dem,4)
     z_fixed = z_L_fixed + z_V_fixed + z_M_fixed + z_A_fixed + z_W_fixed + z_T_fixed;

    z = z_L + z_V + z_M + z_A + z_W + z_T;
    
end
