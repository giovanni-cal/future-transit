%% Single Bus FRT & DRC Comparison
clear
clc
format compact
run Param_def.m

tic %start timer

%LOCAL OPTIMIZATION ROUND
dx = 0.5;
f = 4; % f = 3 corresponds to fixed feeder service, f = 4 to flexible
A = [-1,0,0,1,0,0;0,-1,1,0,0,0];
b = [0;0];
A_eq = []; %A_eq * DX = b_eq
b_eq = [];
lb = [pi/36,0.2,0.1,pi/72,0.03,0.03]; %lower bound of DX
ub = [pi,5,5,pi,1,1]; %upper bound of DX
options = optimoptions('fmincon','Display','off');
options_glob = optimset('Display','off');
DX_mat = []; %matrix with DX(x) in each row;

dr = 0.5;
Z_min = [10^9,10^9,10^9,10^9];
% while dr >= 0.5 & dQ >= 1
     r = 5;
     Z_min(2) = 10^9; %min cost with a given r
     while r < R %optimization procedure varying r
        Q = 0;
        flow_flag = -1;
        Z_min(3) = 10^9; %min cost assuming r and Q, for x > r
        Z_min(4) = 10^9; %min cost assuming r and Q, for x < r
        fval = 0;
        exitflag = 1;
        x = dx/2;
        Z = 0;
        DX_mat = [];
        %Th_r0 = rand*ub(1); S_c0 = rand*ub(2); s0 = rand*ub(3); phi0 = rand*ub(4); H0 = 2*pi/(Th_r0*Q); h0 = 0;
        Th_r0 = (pi/6); S_c0 = 1; s0 = S_c0/2; phi0 = Th_r0/2; H0 = 2*pi/(Th_r0*Q); h0 = H0;
        DX0 = [Th_r0,S_c0,s0,phi0,H0,h0];
        while x < R
            obj_f = @(DX)sum_local_costs(DX,x,r,R,v,v_p,tau_s,tau_p,d,rho_0,gamma,...
            c_L,c_V,c_M,mu_W,mu_A,mu_T,tot_Dem,f);
            f_nlc = @(DX)nl_constr(DX,x,r,Q,R,Cap,tot_Dem,rho_0,gamma,d,f,v,v_p,flow_flag);
            [DX,fval,exitflag] = fmincon(obj_f,DX0,A,b,A_eq,b_eq,lb,ub,f_nlc,options);
            Q = ceil((2*pi)/(DX(1)*DX(5)));
            DX_mat = [DX_mat; x DX Q];
            Z = Z + fval;
            x = x + dx;
        end
        Q_arr = DX_mat(:,end); %array of Q values along x
        Q_0r = Q_arr(1,1);
        Q_rR = min(Q_0r,Q_arr(r/dx+1,1));
        ind1 = 1;
        ind2 = r/dx; %+ ind1% DEMAND(x)
        %capacity constraint on boundary route %tot_Dem*2/pi*(Prob_rR)^2*(H_B/2)/(2*pi) - Cap(1) <= 0
        Prob_rR = integral(Dem_x,r,R)/tot_Dem;
        H_B_up = Cap(1)*4*pi/(tot_Dem*(2/pi)*(Prob_rR)^2);       
    
        while 1
            %local optimization for x > r
            x = dx/2;
            Z = 0;
            Q = Q_0r;
            flow_flag = 1;
            DX_mat = [];
            while x < r %optimal DX values along x
                if x/dx > ind1 %the new Q could decrease towards the periphery
                    flow_flag = 0;
                end
                obj_f = @(DX)sum_local_costs(DX,x,r,R,v,v_p,tau_s,tau_p,d,rho_0,gamma,...
                    c_L,c_V,c_M,mu_W,mu_A,mu_T,tot_Dem,f);
                f_nlc = @(DX)nl_constr(DX,x,r,Q,R,Cap,tot_Dem,rho_0,gamma,d,f,v,v_p,flow_flag);
                [DX,fval,exitflag] = fmincon(obj_f,DX0,A,b,A_eq,b_eq,lb,ub,f_nlc,options);
                if exitflag < 1 %non feasible solution
                    Z = 10^9;
                    break
                end
                Q = (2*pi)/(DX(1)*DX(5));
                DX_mat = [DX_mat; x DX Q];
                Z = Z + fval;
                x = x + dx;
            end
            if exitflag < 1
                if Q_arr((ind1+1),1) >= Q_0r
                    ind1 = ind1 + 1;
                    Q_0r = Q_arr(ind1);
                end
                Q_rR = min(Q_0r,Q_arr(ind2,1));
                continue %skip the remaining part of the while cycle
            end
            %add global costs
            phi_B = DX(4);
            obj_f_glob = @(H_B)sum_global_costs(r,phi_B,H_B,v,tau_s,...
                c_V,c_M,mu_W,mu_A,mu_T,Prob_rR,tot_Dem);
            [H_B,fval,exitflag] = fminbnd(obj_f_glob,0.03,H_B_up,options_glob);
            if exitflag < 1
                continue %skip the remaining part of the while cycle
            end
            Z = Z + fval;
            %local optimization for x > r
            x = r + dx/2;
            %Z = Z_min(4);
            Q = Q_rR;
            flow_flag = 1;
            while x < R %optimal DX values along x
                if x/dx > ind2 %the new Q could decrease towards the periphery
                    flow_flag = 0;
                end
                obj_f = @(DX)sum_local_costs(DX,x,r,R,v,v_p,tau_s,tau_p,d,rho_0,gamma,...
                    c_L,c_V,c_M,mu_W,mu_A,mu_T,tot_Dem,f);
                f_nlc = @(DX)nl_constr(DX,x,r,Q,R,Cap,tot_Dem,rho_0,gamma,d,f,v,v_p,flow_flag);
                [DX,fval,exitflag] = fmincon(obj_f,DX0,A,b,A_eq,b_eq,lb,ub,f_nlc,options);
                if exitflag < 1 %non feasible solution
                    Z = 10^9;
                    break
                end
                Q = (2*pi)/(DX(1)*DX(5));
                DX_mat = [DX_mat; x DX Q];
                Z = Z + fval;
                x = x + dx;
            end
            if Q_arr((ind1+1),1) >= Q_0r
                ind1 = ind1 + 1;
                Q_0r = Q_arr(ind1);
            end
            if ind1 > r/dx
                ind2 = ind2 + 1;
                Q_rR = min(Q_0r,Q_arr(ind2,1));
            end
            if exitflag < 1
                continue %skip the remaining part of the while cycle
            end
            if Z <= Z_min(3) 
                Z_min(3) = Z;
                DX_mat_3 = DX_mat;
                DX_mat = [];
            else
                DX_mat = [];
                break;
            end
        end %optimization with r fixed and Q variable for x > r
        r, DX_mat_3(:,end)
        if Z_min(3) < Z_min(2)
            Z_min(2) = Z_min(3);
            DX_mat_2 = DX_mat_3;
            DX_mat_3 = [];
        else
            DX_mat_3 = [];
            r = r - dr;
            break %exit from the loop and pass to reduced r and Q values
        end
        r = r + dr
     end
%      if Z_min(2) < Z_min(1)
%          Z_min(1) = Z_min(2);
%          DX_mat_1 = DX_mat_2;
%          DX_mat_2 = [];
%      else
%          DX_mat_2 = [];
%          break %exit from the loop and pass to reduced r and Q values
%      end

%end

DX_MAT = DX_mat_2;
DX_MAT(r/dx+1:R/dx,3) = NaN;
DX_MAT(r/dx+1:R/dx,5) = NaN;
DX_MAT(1:r/dx,7) = NaN
%DX_mat(1,1:end)=DX_mat(2,1:end);
Z = 0; Z_L = 0; Z_V = 0; Z_M = 0; Z_A = 0; Z_W = 0; Z_T = 0;
for ii = 1:length(DX_MAT(:,1))
    [z_L,z_V,z_M,z_A,z_W,z_T] = local_cost_fun(DX_MAT(ii,2:end-1),DX_MAT(ii,1),...
        r,R,v,v_p,tau_s,tau_p,d,rho_0,gamma,c_L,c_V,c_M,mu_W,mu_A,mu_T,tot_Dem,f);
     Z_L = Z_L+z_L; Z_V = Z_V+z_V; Z_M = Z_M+z_M;
     Z_A = Z_A+z_A; Z_W = Z_W+z_W; Z_T = Z_T+z_T;
end
[F_V,F_M,F_W,F_T,F_transf] = global_cost_fun(r,phi_B,H_B,v,tau_s,c_V,c_M,...
    mu_W,mu_A,mu_T,Prob_rR,tot_Dem);
 Z_V = Z_V+F_V; Z_M = Z_M+F_M; Z_A = Z_A+F_transf; Z_W = Z_W+F_W; Z_T = Z_T+F_T;
'cost items:'
[Z_L,Z_V,Z_M,Z_A,Z_W,Z_T]
'total cost:'
Z = Z_L + Z_V + Z_M + Z_A + Z_W + Z_T
toc % end timer

figure('DefaultAxesColorOrder', [0 0 1; 1 0 0])
x_ax = linspace(0,R,R/dx);
Th_ax = DX_MAT(1:end,2).*(180/pi);
Sr_ax = DX_MAT(1:end,2).*transpose((dx/2):dx:(R-dx/2));
Sc_ax = DX_MAT(1:end,3);
sr_ax = DX_MAT(1:end,4);
sc_ax = DX_MAT(1:end,5).*transpose((dx/2):dx:(R-dx/2));
H_ax = DX_MAT(1:end,6);
h_ax = DX_MAT(1:end,7);
%subplot(2,2,1);
grid on;
%yticks(0:5);
plot(x_ax,Sr_ax,'-', x_ax,Sc_ax,'-', x_ax,sr_ax,'--', x_ax,sc_ax,'--');
title(['Spacing and Headway | r = ' num2str(r) ' km' ' | R = ' num2str(R) ' km']); % ' | Q = ' num2str(Q) ' veh/h'
xlabel('Distance from center x (km)');
ylabel('Spacing (km)'); %left y-axis
hold on
yyaxis right
plot(x_ax,H_ax,'g',x_ax,h_ax,'g');
ylabel('Headway (hr)'); %right y-axis
legend({'Sr(x)','Sc(x)','sr(x)','sc(x)','H(x)','h(x)'},'Location','northwest','NumColumns',2)
hold off


%ylabel('Spacing between radial lines (km)'); %left y-axis

% 
figure
p = pie([Z_L,Z_V,Z_M,Z_A,Z_W,Z_T]);
labels = {'Z_L','Z_V','Z_M','Z_A','Z_W','Z_T'};
legend(labels,'Location','southoutside','Orientation','horizontal');

if f == 3
    save fixedFMLM.mat
elseif f == 4
    save flexibleFMLM.mat
else
    save onlyBackbone.mat
end
