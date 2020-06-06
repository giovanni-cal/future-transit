function [z] = local_cost_fun(DX,x,r,v,v_p,tau_s,tau_p,d,rho_0,gamma,c_L,c_V,c_M,mu_W,mu_A,mu_T)

Theta_r = DX(1);
S_c     = DX(2);
s       = DX(3);
phi     = DX(4);
H       = DX(5);
    
    %% DEMAND(x)
Dem = 2*pi*x*rho_0.*exp(-gamma*x);
%Prob_x = @(x) Dem/tot_Dem;
%Prob_rR = integral(Prob_x,r,R);
    
    %% COMMERCIAL SPEED
    
%Radial lines
if x > r
    v_comm_r = 1 / (1 / v(1) + (tau_s(1) / s));
else
    v_comm_r = 1 / (1 / v(2) + (tau_s(2) / s));
end

%Ring lines
if x < r
    v_comm_c = 1 / (1 / v(2) + (tau_s(2) / (phi * x)));
end

    % AGENCY COSTS
 %% INFRASTRUCTURE LENGTH
%Local cost for the RAILWAY lines y_L_r(x)
if x > r
    y_L_rail = (2 * pi) / Theta_r;
else
    y_L_rail = 0;
end
%Local cost for the METRO lines y_L_m(x)
if x > r
    y_L_metro = 0;
else
    y_L_metro = (2 * pi) / Theta_r + (2 * pi * x) / S_c;
end
%Local cost for the the feeder BUS lines y_L_b(x)
if x > r
    y_L_bus = ((2 * Theta_r * x - d) / s) * (4 * pi / Theta_r);   %line_length*n_of_lines
else
    y_L_bus = 0;
end

 %% VEHICLE-DISTANCE TRAVELED PER HOUR
%Local cost for the RAILWAY lines y_V_r(x)
if x > r
    y_V_rail = (4 * pi) / (Theta_r * H);
else
    y_V_rail = 0;
end
%Local cost for the METRO lines y_V_m (x)   
if x > r
    y_V_metro = 0;
else
    y_V_metro = (4 * pi) / (Theta_r * H) + (2 * pi * x) / (S_c * H); 
end   
%Local cost for the feeder BUS lines y_V_b (x)   
if x > r
    y_V_bus = (2 * Theta_r * x - d) / H;
else
    y_V_bus = 0;
end
      
 %% FLEET SIZE
%Local cost for the RAILWAY lines y_M_rail(x)
if x > r
    y_M_rail = (4 * pi) / (Theta_r * H * v_comm_r);
else
    y_M_rail = 0;
end 
%Local cost for the METRO lines y_M_metro (x)
if x > r
    y_M_metro = 0;
else
    y_M_metro = (4 * pi) / (Theta_r * H * v_comm_r) + ...
                (4 * pi * x) / (S_c * H * v_comm_c); 
end  
%Local cost for the feeder BUS lines y_M_bus (x)
if x > r
    y_M_bus = ((Theta_r * x - d) / v(3)) + tau_s(3) * (Theta_r * x / d - 1) + ...
              tau_p * Dem * (Theta_r * x - d) / 2 * s * H;
    %         travel time + time lost due to stops + time lost per pax
else
    y_M_bus = 0;
end
    
    %% USER COSTS

%Local cost for the waiting time y_W (x)
if x > r
    y_W = Dem * H;
else
    y_W = Dem * H;
end
%Local cost for the access time y_A (x)
if x > r
    y_A = Dem * (Theta_r*x/4 + s/4)/v_p; %MODIFY!!!
else
    y_A = Dem * (S_c/4 + s/4)/v_p; %MODIFY!!!
end
%Local cost for the in-vehicle time y_T (x)
if x > r
    y_T = Dem * 0.5; %MODIFY!!!
else
    y_T = Dem * 0.5; %MODIFY!!!
end

    %%linear combination
 z_L = c_L(1) * y_L_rail + c_L(2) * y_L_metro + c_L(3) * y_L_bus;
 z_V = c_V(1) * y_V_rail + c_V(2) * y_V_metro + c_V(3) * y_V_bus;
 z_M = c_M(1) * y_M_rail + c_M(2) * y_M_metro + c_M(3) * y_M_bus;
 z_W = mu_W * y_W;
 z_A = mu_A * y_A;
 z_T = mu_T * y_T;

 z = z_L + z_V + z_M + z_W + z_A + z_T;
    
end
