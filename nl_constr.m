function [c,ceq] = nl_constr(DX,x,r,Q,R,Cap,tot_Dem,rho_0,gamma,ns,f,flow_flag)

Th_r = DX(1);
S_c  = DX(2);
s    = DX(3);
H    = DX(5);
h    = DX(6);
d    = DX(7);

Dem = 2*pi*x*rho_0.*exp(-gamma*x); %pax/km
Dem_dy = @(y) (2*pi*y*rho_0.*exp(-gamma*y));
Dem_0x = integral(Dem_dy,0,x)/tot_Dem; %adim
Dem_xR = integral(Dem_dy,x,R)/tot_Dem;


    if x > r
        if f == 1 
            d0 = d/2;%h/(2*(1/v_p-1/v(2)));
            prob_walk = min([(2*d0)/(Th_r*x), 1]);
        elseif f == 2
            d0 = d/2;
            prob_walk = min([(2*d0^2)/(s*Th_r*x), 1]);
        end    
        O_r = tot_Dem*(Dem_xR)*(Th_r*H)/(2*pi); %radial line
        c_r = O_r - Cap(1);
        if f == 1
            O_b = Dem*(s/ns)*Th_r/(4*pi)*h*(1-prob_walk); %FMLM bus
            c_b = O_b - Cap(2);
        elseif f == 2
            O_b = Dem*(s/ns)*Th_r/(4*pi)*h*(1-prob_walk); %FMLM bus
            c_b = O_b - Cap(3);
        else
            c_b = -10000;
        end
        if flow_flag == 1
            c = max([c_r,c_b]);
            ceq = (2*pi)/(Th_r*H) - Q;
        elseif flow_flag == 0
            c = max([c_r,c_b,((2*pi)/(Th_r*H) - Q)]);
            ceq = 0;
        else
            c = max([c_r,c_b]);
            ceq = 0;
        end
    else
        O_r = tot_Dem*(2/pi*Dem_xR*Dem_0x + (1-2/pi)*Dem_xR)*(Th_r*H)/(2*pi); %radial
        c_r = O_r - Cap(1);
        O_c = Dem*2/pi*(Dem_xR)*(S_c*H)/(2*pi); %vehicle occupancy circular line
        c_c = O_c - Cap(1);
        if flow_flag == 1
            c = max([c_r,c_c]);
            ceq = (2*pi)/(Th_r*H) - Q;
        elseif flow_flag == 0
            c = max([c_r,c_c,((2*pi)/(Th_r*H) - Q)]);
            ceq = 0;
        else
            c = max([c_r,c_c]);
            ceq = 0;
        end
    end
    
end