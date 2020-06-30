function [c,ceq] = nl_constr(DX,x,r,Q,R,Cap,tot_Dem,rho_0,gamma,d,f,v,v_p,flow_flag)

Th_r = DX(1);
S_c  = DX(2);
H    = DX(5);
h    = DX(6);

Dem = 2*pi*x*rho_0.*exp(-gamma*x);
Dem_dy = @(y) (2*pi*y*rho_0.*exp(-gamma*y));
Dem_0x = integral(Dem_dy,0,x)/tot_Dem;
Dem_xR = integral(Dem_dy,x,R)/tot_Dem;

if x > r && f < 5
    d1 = (d*v(f)+h/2)/(v(f)-v_p); %fix this value, for each x, for non-peak hours!
    if d1 > (Th_r*x/4) || d1 < d
        d1 = d/2;
    end
end

    if x > r
        O_r = tot_Dem*(Dem_xR)*(Th_r*H)/(2*pi); %radial line
        c_r = O_r - Cap(1);
        if f < 5
            O_b = Dem*(Th_r*h)/(4*pi)*(1-(2*d1-d)/(Th_r*x)); %FMLM bus
            c_b = O_b - Cap(f);
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
        O_c = Dem*S_c*2/pi*(Dem_xR)*(S_c*H)/(2*pi); %vehicle occupancy circular line
        c_c = O_c - Cap(2);
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