%% function that computes the 6 cost components and sums up them in " z ".
% in the ADAPTIVE scenario (f = 3) the function selects the minimum cost
% solution between MRT-only, Fixed or Demand Responsive feeder.
function [z,I] = sum_local_costs_f(DX,x,r,R,v,v_p,tau_s,tau_p,rho_0,gamma,...
    c,mu,tot_Dem,ns,f)

    if f == 0 || f == 1 || f == 2
        [z_L,z_V,z_M,z_A,z_W,z_T] = local_cost_fun(DX,x,r,R,v,v_p,...
            tau_s,tau_p,rho_0,gamma,c,mu,tot_Dem,ns,f);
        z = z_L + z_V + z_M + z_A + z_W + z_T;
        I = f;
        
    elseif f == 3
        
        [z_L_bb,z_V_bb,z_M_bb,z_A_bb,z_W_bb,z_T_bb] = local_cost_fun(DX,x,r,...
            R,v,v_p,tau_s,tau_p,rho_0,gamma,c,mu,tot_Dem,ns,0);
        z_bb = z_L_bb + z_V_bb + z_M_bb + z_A_bb + z_W_bb + z_T_bb;
        
        [z_L_fx,z_V_fx,z_M_fx,z_A_fx,z_W_fx,z_T_fx] = local_cost_fun(DX,x,r,...
            R,v,v_p,tau_s,tau_p,rho_0,gamma,c,mu,tot_Dem,ns,1);
        z_fx = z_L_fx + z_V_fx + z_M_fx + z_A_fx + z_W_fx + z_T_fx;
    
        [z_L_dr,z_V_dr,z_M_dr,z_A_dr,z_W_dr,z_T_dr] = local_cost_fun(DX,x,r,...
            R,v,v_p,tau_s,tau_p,rho_0,gamma,c,mu,tot_Dem,ns,2);
        z_dr = z_L_dr + z_V_dr + z_M_dr + z_A_dr + z_W_dr + z_T_dr;
        
        [z,I] = min([z_bb,z_fx,z_dr]);
        
    else
        error('invalid "f" value');
    end
    
        if x < r
            I = 1;
        end
end