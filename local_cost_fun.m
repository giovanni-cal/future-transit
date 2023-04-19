function [z_L,z_V,z_M,z_A,z_W,z_T] = local_cost_fun(DX,x,r,R,v,v_p,...
    tau_s,tau_p,rho_0,gamma,c,mu,tot_Dem,ns,f)

Th_r = DX(1); %angular spacing btwn radial lines
S_c  = DX(2); %linear spacing btwn circular lines (x < r)
s    = DX(3); %linear stop spacing along radial lines
phi  = DX(4); %angular stop spacing along circular lines (x < r)
H    = DX(5); %MRT headway
h    = DX(6); %feeder headway (x > r)
d    = DX(7); %spacing btwn feeder stops (x > r)

y_L_rail = 0; y_L_metro = 0; y_L_bus = 0; y_L_drt = 0;
y_V_rail = 0; y_V_metro = 0; y_V_bus = 0; y_V_drt = 0;
y_M_rail = 0; y_M_metro = 0; y_M_bus = 0; y_M_drt = 0;
y_A_MRT = 0; y_A_fdr = 0;
y_W_MRT = 0; y_W_fdr = 0;
y_T_MRT = 0; y_T_fdr = 0;

%% DEMAND(x)
Dens = rho_0.*exp(-gamma*x);
Dem_x = Dens * (2*pi*x);
Dem_dy = @(y) (2*pi*y*rho_0.*exp(-gamma*y));
Prob_0x = integral(Dem_dy,0,x)/tot_Dem;
Prob_xR = integral(Dem_dy,x,R)/tot_Dem;

%% COMMERCIAL SPEED

%Radial lines
v_comm_r = 1/(1/v(1)+(tau_s(1)/s));

%Ring lines
v_comm_c = 1/(1/v(1)+(tau_s(1)/(phi*x)));

%% FEEDER BUS CYCLE AND WALK PROBABILITY

tau_T = 0.025;  %idle-time at terminal
if x > r
    if ns == 1
        dL = 0;
    else
        dL = s/4; %additional distance due to the strips
    end
    if f == 0
        prob_walk = 1; %probability that a user walk directly to the station
    elseif f == 1
        d0 = d/2;%h/(2*(1/v_p-1/v(2)));
        prob_walk = min([(2*d0)/(Th_r*x), 1]); %ratio btwn the walking area and the service region
        nn = Dens*(1-prob_walk)*h*(s/ns)*Th_r*x; %it is 2*Dens*...*Th_r*x/2
        if ns > 1
            CL_fx = (s/2 + Th_r*x-d);
        else
            CL_fx = (Th_r*x-d);
        end
        C_fx = tau_T + (CL_fx/v(2)) + ... % terminal time + travel time
            tau_s(2)*(Th_r*x/d - 1) + tau_p*nn; % + time lost per stop + per pax
        C_fx = max(C_fx,h);
    elseif f == 2
        d0 = d/2; %min([h/(2*(1/v_p-1/v(3))),s/2,Th_r*x/2]);
        %s < 0.5 km --> d0 < 0.25 km
        prob_walk = min([(2*d0^2)/(s*Th_r*x), 1]);
        nn = Dens*(1-prob_walk)*h*(s/ns)*Th_r*x; %it is 2*Dens*...*Th_r*x/2 
        %CL = Th_r*x*(nn/(nn+1)) + s/6*nn + 2/3*s; %Quadrifoglio&Li (2009)
        CL_dr = Th_r*x*(nn/(nn+1)) + (s/ns)/3*nn + (s/ns)/2; %Quadrifoglio (2006)
        C_dr = tau_T + CL_dr/v(3) + nn * (tau_s(3)+tau_p);
        C_dr = max(C_dr,h);
    else
        error('invalid "f" value');
    end
end

%% AGENCY LOCAL COSTS

% INFRASTRUCTURE LENGTH (km/km -> adim.)
% Capital cost 

if x > r
    y_L_rail = 2*pi/Th_r;
    if f == 0
        %do nothing
    elseif f == 1
        y_L_bus = ns * (Th_r*x)/s * (2*pi/Th_r);   %line_length*n_of_lines)
    elseif f == 2
        %the same of y_L_bus since no new infrastructure is needed for DRT
        y_L_drt = ns * (Th_r*x)/s * (2*pi/Th_r);
    else
        error('invalid "f" value');
    end
    c_st = 100/s;
else
    y_L_metro = 2*pi/Th_r + 2*pi*x/S_c;
    c_st = 300/s + 300/(phi*x);
end

% VEHICLE-DISTANCE TRAVELED PER HOUR ((veh*km/h)/km -> veh/h)
% Operation cost

if x > r
    y_V_rail = (4*pi)/(Th_r*H);
    if f == 0
        %do nothing
    elseif f == 1
        y_V_bus = ns * CL_fx/(s*h) * (4*pi/Th_r);
    elseif f == 2
        y_V_drt = ns * CL_dr/(s*h) * (4*pi/Th_r);
    else
        error('invalid "f" value');
    end
else
    y_V_metro = (4*pi)/(Th_r*H) + (2*pi*x)/(S_c*H);
end

% FLEET SIZE ((veh)/km -> veh/km)
% Vehicles & crew cost

if x > r
    y_M_rail = (4*pi)/(Th_r*H*v_comm_r);
    if f == 0
        %do nothing
    elseif f == 1
        y_M_bus = ns * C_fx/(s*h) * (4*pi/Th_r);
    elseif f == 2
        y_M_drt = ns * C_dr/(s*h) * (4*pi/Th_r);
    else
        error('invalid "f" value');
    end
else
    y_M_metro = (4*pi)/(Th_r*H*v_comm_r) + (4*pi*x)/(S_c*H*v_comm_c);
end

%% USER LOCAL COSTS 

%Local cost for the access time y_A (x) (pax*h/km)

if x > r
    if f == 0
        y_A_MRT = 2*Dem_x*(Th_r*x/4+s/4)/v_p;
    elseif f == 1
        y_A_MRT = 2*Dem_x*(1/v_p*(d0/2+s/4)*(prob_walk));
        y_A_fdr = 2*Dem_x*(1/v_p*(d/4+(s/ns)/4)*(1-prob_walk)); %2*Dem_x*1/v_p*((Th_r/4+s/4));
    elseif f == 2
        %if d0 > s/2
            %avg_walk_dist = (2*d0^2-d0*s+s^2/6)/(4*d0-s) + (3*d0*s-s^2)/(12*d0-3*s);
        %else
        avg_walk_dist = d0/3 + d0/3;
        %end
    	y_A_MRT = 2*Dem_x*(prob_walk)*avg_walk_dist/v_p;
    else
        error('invalid "f" value');
    end
else
    y_A_MRT = 2*Dem_x *(Prob_0x*(Th_r*x/4 + s/4)/v_p ...
        + Prob_xR*(S_c/4 + phi*x/4)/v_p);
end

%Local cost for the waiting time y_W (x) (pax*h/km)

W_red_0R = Prob_0x*Th_r/pi;
W_red_rR = (Dem_x*s/tot_Dem)*Th_r/pi;
if x > r
    if f == 0
        %2 times Dem_x for an avg wt time of H/2
        y_W_MRT = Dem_x * H * (1-W_red_0R-W_red_rR);
    elseif f == 1
        y_W_MRT = Dem_x * H * (1-W_red_0R-W_red_rR);
        y_W_fdr = Dem_x * h * (1-prob_walk); %* (1-W_red_0R-W_red_rR);
    elseif f == 2
        y_W_MRT = Dem_x * H * (1-W_red_0R-W_red_rR);
        y_W_fdr = Dem_x * h * (1-prob_walk); %* (1-W_red_0R-W_red_rR);
    else
        error('invalid "f" value');
    end
else
    %2 times H/2 for each user
    y_W_MRT = Dem_x * H * (1-W_red_0R); 
end

%Local cost for the in-vehicle time y_T (x) (pax*h/km)

T_red_rR = (Prob_0x*Prob_xR + (Dem_x*s/tot_Dem)^2) * (Th_r/pi);
if x > r
    if f == 0
        y_T_MRT = tot_Dem*2*(Prob_xR-T_red_rR)/v_comm_r;
    elseif f == 1
        y_T_MRT = tot_Dem*2*(Prob_xR-T_red_rR)/v_comm_r;
        y_T_fdr = 2*Dem_x*(C_fx/4+(d0/2+dL)/v(2))*(1-prob_walk);
    elseif f == 2
    	%y_T = tot_Dem*2*Prob_xR*1/v_comm_r + 2*Dem_x*(C/2)*(1-prob_walk);
        y_T_MRT = tot_Dem*2*(Prob_xR-T_red_rR)/v_comm_r;
        y_T_fdr = 2*Dem_x*(C_dr/4+dL/v(3))*(1-prob_walk);
    else
        error('invalid "f" value');
    end
else
    y_T_MRT = 2*Dem_x*Prob_xR*(2/pi)*(x/v_comm_c)...
        + tot_Dem*(2*Prob_xR*(1-2/pi)+2*Prob_0x*Prob_xR*(2/pi))/v_comm_r;
end

%%linear combination
z_L = (0.5*c(1,1)+c_st)*y_L_rail + (c(1,1)+c_st)*y_L_metro + c(1,2)*y_L_bus + c(1,3)*y_L_drt;
z_V = c(2,1)*y_V_rail + c(2,1)*y_V_metro + c(2,2)*y_V_bus + c(2,3)*y_V_drt;
z_M = c(3,1)*y_M_rail + c(3,1)*y_M_metro + c(3,2)*y_M_bus + c(3,3)*y_M_drt;
z_A = mu(3,1) * (y_A_MRT + y_A_fdr);
z_W = mu(2,1) * (y_W_MRT + y_W_fdr);
z_T = mu(1,1) * (y_T_MRT + y_T_fdr);
 
end
