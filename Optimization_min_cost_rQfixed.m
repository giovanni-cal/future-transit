clear
clc
format compact
run Param_def.m

tic %start timer

%LOCAL OPTIMIZATION ROUND
dx = 0.5;
Q = 90;
r = 5;
f = 3; %remove
A = [-1,0,0,1,0,0;0,-1,1,0,0,0];
b = [0;0];
A_eq = []; %A_eq * DX = b_eq
b_eq = [];
lb = [pi/36,0.2,0.1,pi/72,0.01,0.01]; %lower bound of DX
ub = [pi,5,5,pi,1,1]; %upper bound of DX
DX_mat = []; %matrix with DX(x) in each row;

%Th_r0 = rand*ub(1); S_c0 = rand*ub(2); s0 = rand*ub(3); phi0 = rand*ub(4); H0 = 2*pi/(Th_r0*Q); h0 = H0;
Th_r0 = (pi/6); S_c0 = 1; s0 = S_c0/2; phi0 = Th_r0/2; H0 = 2*pi/(Th_r0*Q); h0 = H0;
DX0 = [Th_r0,S_c0,s0,phi0,H0,h0];
x = dx/2;
while x < R %optimal DX values along x
    obj_f = @(DX)sum_local_costs(DX,x,r,R,v,v_p,tau_s,tau_p,d,rho_0,gamma,...
        c_L,c_V,c_M,mu_W,mu_A,mu_T,tot_Dem,f);
    f_nlc = @(DX)nl_constr(DX,x,r,Q,R,Cap,tot_Dem,rho_0,gamma,d,f,v,v_p);

    [DX,fval,exitflag] = fmincon(obj_f,DX0,A,b,A_eq,b_eq,lb,ub,f_nlc);

    if exitflag < 1 %non feasible solution
        break
    end

    if x < r
        DX(6) = 0;
    else
        DX(2) = 0; DX(4) = 0;
    end

    %Q = (2*pi)/(DX(1)*DX(5));
    DX_mat = [DX_mat; x DX (2*pi)/(DX(1)*DX(5))];
    x = x + dx;
end %optimization with r and Q fixed

DX_mat
%DX_mat(1,1:end)=DX_mat(2,1:end);
Z = 0; Z_L = 0; Z_V = 0; Z_M = 0; Z_A = 0; Z_W = 0; Z_T = 0;
for ii = 1:1:length(DX_mat(:,1))
    [z_L,z_V,z_M,z_A,z_W,z_T] = local_cost_fun(DX_mat(ii,2:end),DX_mat(ii,1),r,R,v,v_p,tau_s,tau_p,d,...
        rho_0,gamma,c_L,c_V,c_M,mu_W,mu_A,mu_T,tot_Dem,f);
     Z_L = Z_L+z_L; Z_V = Z_V+z_V; Z_M = Z_M+z_M;
     Z_A = Z_A+z_A; Z_W = Z_W+z_W; Z_T = Z_T+z_T;
     Z = Z + z_L + z_V + z_M + z_A + z_W + z_T;
end
'cost items:'
[Z_L,Z_V,Z_M,Z_A,Z_W,Z_T]
'total cost:'
Z
toc % end timer

figure
lin_sp = DX_mat(1:end,2).*transpose((dx/2):dx:(R-dx/2));
[plotTH,Th_r_line,H_line] = plotyy((dx/2):dx:(R-dx/2),lin_sp,...
    (dx/2):dx:(R-dx/2),DX_mat(1:end,6));
yticks(0:5);
grid on;
% axis(plotTH(1),[0 R 0 10])
% axis(plotTH(2),[0 R 0 max(DX_mat(:,5))+0.1])
title(['Spacing and Headway | r = ' num2str(r) ' km' ' | R = ' num2str(R) ' km' ' | Q = ' num2str(Q) ' veh/h']);
xlabel('Distance from center x (km)');
ylabel(plotTH(1),'Spacing between radial lines (km)'); %left y-axis
ylabel(plotTH(2),'Headway (hr)'); %right y-axis
% 
figure
p = pie([Z_L,Z_V,Z_M,Z_A,Z_W,Z_T]);
labels = {'Z_L','Z_V','Z_M','Z_A','Z_W','Z_T'};
legend(labels,'Location','southoutside','Orientation','horizontal');

