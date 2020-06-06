function [c,ceq] = nl_constr(DX,Q,tot_Dem,Cap) %x

Theta_r = DX(1);
H       = DX(5);

    c = tot_Dem * 0.001 * (Theta_r * H) / (2 * pi) - Cap;
    ceq = (2 * pi) / (Theta_r * H) - Q;
    
end