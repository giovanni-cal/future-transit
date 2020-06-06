%% Single Bus FRT & DRC Comparison
clear
clc
format compact
run Param_def.m

tic %start timer

%LOCAL OPTIMIZATION ROUND
dx = 0.25;
r = 8;
Q = 10;
R = 25; %remove
lb = [0,0.1,0.1,0,0.033];
ub = [pi,r,r,pi,1];
DX_mat = [0,0,0,0,0]; %matrix with DX(x) in each row
DX0 = [(pi/6),1,0.5,(pi/12),0.1]; %DX0 = [Theta_r,S_c,s,phi,H]
for x = (dx/2):dx:(R-dx/2)
    obj_f = @(DX)local_cost_fun(DX,x,r,v,v_p,tau_s,tau_p,d,rho_0,gamma,c_L,c_V,c_M,mu_W,mu_A,mu_T);
    f_nlc = @(DX)nl_constr(DX,Q,tot_Dem,Cap);
    
    DX = fmincon(obj_f,DX0,[],[],[],[],lb,ub,f_nlc);
    DX_mat = [DX_mat;DX];
end

DX_mat(1,1:end)=DX_mat(2,1:end);
DX_mat
figure
[plotTH,Theta_r_line,H_line] = plotyy(0:dx:R,DX_mat(1:end,1),0:dx:R,DX_mat(1:end,5));
axis(plotTH(1),[0 R 0 4])
axis(plotTH(2),[0 R 0 1])
title('Spacing and Headway');
xlabel('Distance from center x');
ylabel(plotTH(1),'Theta_r'); %left y-axis
ylabel(plotTH(2),'H'); %right y-axis

toc % end timer