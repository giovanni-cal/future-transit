function [z_L,z_V,z_M,z_A,z_W,z_T] = local_cost_fun(DX,x,r,R,v,v_p,...
    tau_s,tau_p,d,rho_0,gamma,c_L,c_V,c_M,mu_W,mu_A,mu_T,tot_Dem,f)

Th_r = DX(1);
S_c  = DX(2);
s    = DX(3);
phi  = DX(4);
H    = DX(5);
h    = DX(6);

%% DEMAND(x)
Dem = 2*pi*x*rho_0.*exp(-gamma*x);
Dem_dy = @(y) (2*pi*y*rho_0.*exp(-gamma*y));
Dem_0x = integral(Dem_dy,0,x)/tot_Dem;
Dem_xR = integral(Dem_dy,x,R)/tot_Dem;

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

%% FEEDER BUS CYCLE
if x > r
    if f == 3
        C = ((Th_r*x - d)/v(f)) + tau_s(f)*(Th_r*x/d - 1) + tau_p*Dem*(Th_r*x - d)/2*s*H;
        %       travel time     + time lost due to stops  + time lost per pax
    elseif f == 4
        nn = 2*Dem*h*s*Th_r*x/2;
        CL = Th_r*x/2*(1-1/(nn+1)) + 2*(nn-1)*s/3 + s/2;
        C = CL/v(f) + nn * (tau_s(f)+tau_p);
    else
        C = 0;
    end
end

% AGENCY COSTS

%% INFRASTRUCTURE LENGTH
%Local cost for the RAILWAY lines y_L_r(x)
if x > r
    y_L_rail = 2*pi/Th_r;
else
    y_L_rail = 0;
end
%Local cost for the METRO lines y_L_m(x)
if x > r
    y_L_metro = 0;
else
    y_L_metro = 2*pi/Th_r + 2*pi*x/S_c;
end
%Local cost for the the feeder BUS lines y_L_b(x)
% if x > r
%     y_L_bus = ((Th_r*x - d)/s) * (4*pi/Th_r);   %line_length*n_of_lines
% else
%     y_L_bus = 0;
% end

%% VEHICLE-DISTANCE TRAVELED PER HOUR
%Local cost for the RAILWAY lines y_V_r(x)
if x > r
    y_V_rail = (4*pi)/(Th_r*H);
else
    y_V_rail = 0;
end
%Local cost for the METRO lines y_V_m (x)
if x > r
    y_V_metro = 0;
else
    y_V_metro = (4*pi)/(Th_r*H) + (2*pi*x)/(S_c*H);
end
%Local cost for the feeder BUS lines y_V_b (x)
if x > r
    if f == 3
        y_V_bus = (Th_r*x - d)/h*(4*pi/Th_r);
    elseif f == 4
        y_V_bus = CL/h * (4*pi/Th_r);
    else
        y_V_bus = 0;
    end
else
    y_V_bus = 0;
end

%% FLEET SIZE
%Local cost for the RAILWAY lines y_M_rail(x)
if x > r
    y_M_rail = (4*pi)/(Th_r*H*v_comm_r);
else
    y_M_rail = 0;
end
%Local cost for the METRO lines y_M_metro (x)
if x > r
    y_M_metro = 0;
else
    y_M_metro = (4*pi)/(Th_r*H*v_comm_r) + (4*pi*x)/(S_c*H*v_comm_c);
end
%Local cost for the feeder BUS lines y_M_bus (x)
if x > r
    y_M_bus = (C/h)*(4*pi/Th_r);
else
    y_M_bus = 0;
end

%% USER COSTS

if x > r
    if f < 5
        d1 = (d*v(f)+h/2)/(v(f)-v_p); %fix this value, for each x, for non-peak hours!
        if d1 > (Th_r*x/4) || d1 < d
            d1 = d/2;
        end
        prob_walk = (2*d1-d)/(Th_r*x); %probability that a user walk directly to the station instead of taking the bus
    else
        prob_walk = 1; %probability that a user walk directly to the station instead of taking the bus
    end
end
 

%Local cost for the access time y_A (x)
if x > r
    if f == 3
        y_A = 2*Dem*1/v_p*(((2*d1-d)/4+s/4)*(prob_walk) +...
            (d/4+s/4)*(1-prob_walk)); %2*Dem*1/v_p*((Th_r/4+s/4));
    elseif f == 4
    	y_A = 2*Dem*1/v_p*((2*d1-d)/4+s/4)*(prob_walk); %modify!
    else
        y_A = 2*Dem*1/v_p*(Th_r*x/4+s/4);
    end
else
    y_A = 2*Dem*1/v_p*(2/pi*(1/4*(Th_r*x + s)*Dem_0x + 1/4*(S_c + phi*x)*Dem_xR)...
        + (1-2/pi) * 1/4*(Th_r*x + s));
end
%Local cost for the waiting time y_W (x)
if x > r
    if f == 3
        y_W = Dem * (H + h*(1-prob_walk));
    elseif f == 4
        y_W = Dem * (H + h*(1-prob_walk)); %modify!
    else
        y_W = Dem * H;
    end
else
    y_W = Dem * H; %2 times H/2 for each user
end
%Local cost for the in-vehicle time y_T (x)
if x > r
    if f == 3
        y_T = tot_Dem*2*Dem_xR*1/v_comm_r + 2*Dem*(C/4+(d1-d/2)/(2*v(f)))*(1-prob_walk);
    elseif f == 4
    	y_T = tot_Dem*2*Dem_xR*1/v_comm_r + 2*Dem*(C/4+(d1-d/2)/(2*v(f)))*(1-prob_walk); %modify!
    else
        y_T = tot_Dem*2*Dem_xR*1/v_comm_r;
    end
else
    y_T = tot_Dem*2*((1-2/pi)+2/pi*Dem_0x)*Dem_xR*1/v_comm_r...
        + Dem*2*2/pi*Dem_xR*x/v_comm_c;
end

%%linear combination
z_L = c_L(1) * y_L_rail + c_L(2) * y_L_metro; %+ c_L(f) * y_L_bus;
z_V = c_V(1) * y_V_rail + c_V(2) * y_V_metro + c_V(f) * y_V_bus;
z_M = c_M(1) * y_M_rail + c_M(2) * y_M_metro + c_M(f) * y_M_bus;
z_A = mu_A * y_A;
z_W = mu_W * y_W;
z_T = mu_T * y_T;

end
