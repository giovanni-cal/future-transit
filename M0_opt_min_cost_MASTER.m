%master script that optimizes the transit network structure
clc
clear
format compact
%doplots = 1;

tic %start timer

%myinput can be taken from a txt file ore set manually
%read in this order: R ; rho_0 ; gamma ; f ; peak_flag.
%myinput = fscanf(fopen('myinput.txt'),'%f',[4,Inf]); fclose('all');
myinput = [25;1600;0.12;1];

% f = 0 -> MRT Only, =1 -> MRT-FRF scheme, =2 -> MRT-DRF scheme, =3 -> adaptive scheme
for n_run = 1:length(myinput(1,:))

%the following parameters can be taken from a txt file or set manually
    R = myinput(1,n_run);
    rho_0 = myinput(2,n_run);
    gamma = myinput(3,n_run);
    f = myinput(4,n_run);
    run M1_Param_def.m %procedure to set the input parameters
    
    DX0 = lb; 
%     Th_r0 = lb(1)+rand*(ub(1)-lb(1));
%     S_c0 = lb(2)+rand*(ub(2)-lb(2));
%     s0 = lb(3)+rand*(ub(3)-lb(3));
%     phi0 = lb(4)+rand*(ub(4)-lb(4));
%     H0 = lb(5)+rand*(ub(5)-lb(5));
%     h0 = lb(6)+rand*(ub(6)-lb(6));
%     DX0 = [Th_r0,S_c0,s0,phi0,H0,h0];
    dx = 1; %%%
    dr = 1; %%%
    r_start = 5%ceil(R/5)
    r_opt = r_start;
    Z_min = 10^9.*ones(2,1); %%%%%

    while dr >= 0.5 %%%%%
        DX_mat_2 = zeros(R/dx,10);
        r = r_start;
        Z_min(1,:) = 10^9; %min cost with a given r
        %optimization procedure varying r
        while r < R
            Q = 0;
            Z_min(2,:) = 10^9; %min cost assuming r and Q, for x > r
            fval = 0;
            exitflag = 1;
            x = dx/2;
            Z = 0;
            Z_arr = zeros(R/dx,1);
            DX_mat = zeros(R/dx,10); %matrix that will be filled with DX(x) in each row;
            countdx = 1;
            DX0var = DX0;
            %capacity constraint on boundary route %tot_Dem*2/pi*(Prob_rR)^2*(H_B/2)/(2*pi) - Cap(1) <= 0
            Prob_rR = integral(Dem_x,r,R)/tot_Dem;
            H_B_up = Cap(1)*4*pi/(tot_Dem*(2/pi)*(Prob_rR)^2);       
            
            Q_try = 450 %%%
            Q_try_best = Q_try;
            Q_step = 100; %%%
            dQ = 1;

            while 1
                x = dx/2;
                Z = 0;
                Z_arr = zeros(R/dx,1);
                Q = Q_try;
                flow_flag = 1;
                DX_mat = zeros(R/dx,10);
                countdx = 1;
                DX0var = DX0;
                lbvar = lb;
                while x < r %optimal DX values along x
                    ns = 1;
                    obj_f = @(DX)sum_local_costs(DX,x,r,R,v,v_p,...
                        tau_s,tau_p,rho_0,gamma,c,mu,tot_Dem,ns,0);
                    f_nlc = @(DX)nl_constr(DX,x,r,Q,R,Cap,tot_Dem,rho_0,...
                        gamma,ns,0,flow_flag);
                    [DX,fval,exitflag] = fmincon(obj_f,DX0var,A,b,A_eq,b_eq,...
                        lbvar,ub,f_nlc,options);
                    if exitflag < 1 %non feasible solution
                        Z = 10^9;
                        break
                    end
                    Q = round((2*pi)/(DX(1)*DX(5)));
                    DX_mat(countdx,:) = [x DX ns Q];
                    Z = Z + dx*fval;
                    Z_arr(countdx,1) = dx*fval;
                    x = x + dx;
                    countdx = countdx + 1;
                    DX0var = DX0;
                    lbvar(5) = DX(5);
                end
                if exitflag < 1
                    Q_try = Q_try + dQ*Q_step
                    continue %skip the remaining part of the while cycle
                end
                phi_B = DX(4);
                %local optimization for x > r
                x = r + dx/2;
                Q = Q_try;
                ns_max = 4;
                flow_flag = 1;
                %%% start of the variables optimization in the suburbs %%%
                while x < R %optimal DX values along x
                    if x > r %the new Q could decrease towards the periphery
                        flow_flag = 0;
                    end
                    ns = 1;
                    if f == 0
                        obj_f = @(DX)sum_local_costs(DX,x,r,R,v,v_p,tau_s,tau_p,...
                            rho_0,gamma,c,mu,tot_Dem,ns,f);
                        f_nlc = @(DX)nl_constr(DX,x,r,Q,R,Cap,tot_Dem,...
                            rho_0,gamma,ns,f,flow_flag);
                        [DX,fval,exitflag] = fmincon...
                            (obj_f,DX0var,A,b,A_eq,b_eq,lbvar,ub,f_nlc,options);
                    elseif f == 1 || f == 2
                        fval_try = 10^9;
                        fval = fval_try;
                        DX_try = DX0;
                        while fval_try <= fval
                            fval = fval_try;
                            DX = DX_try;
                            obj_f = @(DX)sum_local_costs(DX,x,r,R,v,v_p,tau_s,tau_p,...
                                rho_0,gamma,c,mu,tot_Dem,ns,f);
                            f_nlc = @(DX)nl_constr(DX,x,r,Q,R,Cap,tot_Dem,...
                                rho_0,gamma,ns,f,flow_flag);
                            [DX_try,fval_try,exitflag_try] = fmincon...
                                (obj_f,DX0var,A,b,A_eq,b_eq,lbvar,ub,f_nlc,options);
                            if exitflag_try < 1
                                fval_try = 10^9;
                            end
                            if fval_try > fval || ns > ns_max
                                ns = ns -1;
                                break
                            end
                            ns = ns + 1;
                        end
                    elseif f == 3 %adaptive structure
                        obj_f = @(DX)sum_local_costs(DX,x,r,R,v,v_p,tau_s,tau_p,...
                            rho_0,gamma,c,mu,tot_Dem,ns,0);
                        f_nlc = @(DX)nl_constr(DX,x,r,Q,R,Cap,tot_Dem,...
                            rho_0,gamma,ns,0,flow_flag);
                        [DX_bb,fval_bb,exitflag_bb] = fmincon...
                            (obj_f,DX0var,A,b,A_eq,b_eq,lbvar,ub,f_nlc,options);
                        if exitflag_bb < 1
                            fval_bb = 10^9;
                        end
                        %optimizing the number of strips ns_fx
                        ns = 1;
                        fval_fx_try = 10^9;
                        fval_fx = fval_fx_try;
                        DX_fx_try = DX0;
                        while fval_fx_try <= fval_fx
                            fval_fx = fval_fx_try;
                            DX_fx = DX_fx_try;
                            obj_f = @(DX)sum_local_costs(DX,x,r,R,v,v_p,tau_s,tau_p,...
                                rho_0,gamma,c,mu,tot_Dem,ns,1);
                            f_nlc = @(DX)nl_constr(DX,x,r,Q,R,Cap,tot_Dem,...
                                rho_0,gamma,ns,1,flow_flag);
                            [DX_fx_try,fval_fx_try,exitflag_fx_try] = fmincon...
                                (obj_f,DX0var,A,b,A_eq,b_eq,lbvar,ub,f_nlc,options);
                            if exitflag_fx_try < 1
                                fval_fx_try = 10^9;
                            end
                            if fval_fx_try > fval_fx || ns > ns_max
                                ns = ns - 1;
                                break
                            end
                            ns = ns + 1;
                        end
                        %choosing the optimal FMLM setting
                        [fval,index] = min([fval_bb,fval_fx]);
                        if index == 1
                            DX = DX_bb;
                        elseif index == 2
                            DX = DX_fx;
                        else
                            error('invalid DX');
                        end          
                    else
                        error('invalid "f" value');
                    end     
                    if exitflag < 1 || fval >= 10^9 %non feasible solution
                        Z = 10^9;
                        break
                    end
                    Q = round((2*pi)/(DX(1)*DX(5)));
                    DX_mat(countdx,:) = [x DX ns Q];
                    ns_max = min([floor(DX(1,3)/Wstr_max),ns]); %%%
                    Z = Z + dx*fval;
                    Z_arr(countdx,1) = dx*fval;
                    x = x + dx;
                    countdx = countdx + 1;
                    DX0var = DX;
                    lbvar(5) = DX(5);
                end
                if f == 3 && Z < 10^9
                    x = x - dx;
                    countdx = countdx - 1;
                    while x > r
                        Q = DX_mat(countdx,end);
                        ns = DX_mat(countdx,end-1);
                        DX0var = DX_mat(countdx,2:8);
                        lbvar(5) = DX_mat(countdx,6); %the H for the fixed scenario
                        fval_try = 10^9;
                        fval = fval_try;
                        DX_try = DX0;
                        while fval_try <= fval
                            fval = fval_try;
                            DX = DX_try;
                            obj_f = @(DX)sum_local_costs(DX,x,r,R,v,v_p,tau_s,tau_p,...
                                rho_0,gamma,c,mu,tot_Dem,ns,2);
                            f_nlc = @(DX)nl_constr(DX,x,r,Q,R,Cap,tot_Dem,...
                                rho_0,gamma,ns,2,flow_flag);
                            [DX_try,fval_try,exitflag_try] = fmincon...
                                (obj_f,DX0var,A,b,A_eq,b_eq,lbvar,ub,f_nlc,options);
                            if exitflag_try < 1
                                fval_try = 10^9;
                            end
                            ns_max = (DX_try(1,3)/Wstr_max); %%%
                            if fval_try > fval || ns > ns_max
                                ns = ns - 1;
                                break
                            end
                            ns = ns + 1;
                        end
                        if fval < Z_arr(countdx,1)/dx
                            Z = Z + dx*fval - Z_arr(countdx,1);
                            Z_arr(countdx,1) = dx*fval;
                            DX_mat(countdx,:) = [x DX ns Q];
                        else
                            break
                        end
                        x = x - dx;
                        countdx = countdx - 1;
                    end
                    fval = 0;
                end
                if fval >= 10^9
                    Q_try = Q_try + dQ*Q_step
                    continue %skip the remaining part of the while cycle
                end
                %%% end of the varible optimization in the suburbs %%%
                Th_r_min = min(DX_mat(:,2));
                %add global costs
                obj_f_glob = @(H_B)sum_global_costs(r,phi_B,H_B,Th_r_min,...
                    v,tau_s,c,mu,Prob_rR,tot_Dem);
                [H_B,fval,exitflag] = fminbnd(obj_f_glob,0.042,H_B_up,options_glob);
                if exitflag < 1
                    error('global cost function does not satisfy constraints');
                end
                Z = Z + fval;
                if exitflag < 1
                    Q_try = Q_try + dQ*Q_step
                    continue %skip the remaining part of the while cycle
                end
                %if the total cost increases, the optimal Q is found.
                if Z > mean(Z_min(2,:))
                    Q_try = Q_try_best - dQ*Q_step/2
                    switch dQ
                        case 1
                            dQ = 0.2;
                        case 0.2
                            dQ = 0.05;
        %                  case 0.05
        %                      dQ = 0.025;
                        otherwise
                            break; %end of iterations changing Q_try
                    end
                    Z_min(2,:) = 10^9;
                    continue;
                end
                if Z <= min(Z_min(2,:))
                    DX_mat_2 = DX_mat; %new min total cost. Update decision variables matrix.
                    Q_try_best = Q_try
                    H_B_opt = H_B;
                end
                for jj = 1:length(Z_min(2,:))-1
                    Z_min(2,jj) = Z_min(2,jj+1);
                end
                Z_min(2,end) = Z
                disp('********************')
                Q_try = Q_try + dQ*Q_step
            end %optimization with r fixed and Q variable for x > r
%             if min(Z_min(2,:)) > mean(Z_min(1,:))
%                 r_opt = r - dr; %the previous r was the best 
%                 break %exit the loop and reduce dr and dx for a more accurate analysis
%             end
            %end of 2nd iteration with dr = 0.5 km
            if min(Z_min(2,:)) <= min(Z_min(1,:))
                DX_mat_1 = DX_mat_2;
                r_opt = r;
                H_B_best = H_B_opt;
                break %exit the loop and reduce dr and dx for a more accurate analysis
            end
            for jj = 1:length(Z_min(1,:))-1
                Z_min(1,jj) = Z_min(1,jj+1);
            end
            Z_min(1,end) = min(Z_min(2,:));
            disp('##############################')
            if dr <= 0.5 && r == r_start + 2*dr
                break
            end
            r = r + dr %try a higher r
        end
        dr = dr/2;
        [Z_min_fin,pos] = min(Z_min(1,:));
        if dr >= 0.5 %%%%%
            dx = dx/2;
            disp('##############################')
            r_start  = r+2*(pos-length(Z_min(1,:)))*dr - 3*dr
        end
    end
    DX_mat = DX_mat_1;
    DX_mat(r_opt/dx+1:R/dx,3) = NaN;
    DX_mat(r_opt/dx+1:R/dx,5) = NaN;
    DX_mat(1:r_opt/dx,7) = NaN;
    DX_mat(1:r_opt/dx,8) = NaN;
    if f == 0
        DX_mat(:,7) = NaN;
        DX_mat(:,8) = NaN;
    end

    switch f
        case 0
            FMLM = [DX_mat(:,1),zeros(length(DX_mat(:,1)),1)];
        case 1
            FMLM = [DX_mat(:,1),ones(length(DX_mat(:,1)),1)];
        case 2
            FMLM = [DX_mat(:,1),1+ones(length(DX_mat(:,1)),1)];
        case 3
            FMLM = [DX_mat(:,1),zeros(length(DX_mat(:,1)),1)];
            for ii = 1:length(FMLM)
                [z,I] = sum_local_costs_f(DX_mat(ii,2:end-2),DX_mat(ii,1),r_opt,R,...
                    v,v_p,tau_s,tau_p,rho_0,gamma,c,mu,tot_Dem,DX_mat(ii,end-1),f);
                FMLM(ii,2) = I-1;
            end
        otherwise
            error('invalid "f" value');
    end

    Z_mat = zeros(length(DX_mat(:,1)),8);
    for ii = 1:length(DX_mat(:,1))
        Z_mat(ii,1) = ii*dx-dx/2;
        [Z_mat(ii,2),Z_mat(ii,3),Z_mat(ii,4),Z_mat(ii,5),Z_mat(ii,6),Z_mat(ii,7)]...
            = local_cost_fun(DX_mat(ii,2:end-2),DX_mat(ii,1),...
            r_opt,R,v,v_p,tau_s,tau_p,rho_0,gamma,c,mu,tot_Dem,DX_mat(ii,end-1),FMLM(ii,2));
        Z_mat(ii,2:7) = dx * Z_mat(ii,2:7);
        Z_mat(ii,end) = integral(Dem_x,Z_mat(ii,1)-dx/2,Z_mat(ii,1)+dx/2);
    end
    phi_B = DX_mat(r_opt/dx,5);
    Th_r_min = min(DX_mat(:,2));
    H_B = H_B_best;
    [F_V,F_M,F_W,F_T,F_transf] = ...
        global_cost_fun(r_opt,phi_B,H_B_best,Th_r_min,v,tau_s,c,mu,Prob_rR,tot_Dem);
    Z_L = sum(Z_mat(:,2));
    Z_V = sum(Z_mat(:,3))+F_V;
    Z_M = sum(Z_mat(:,4))+F_M;
    Z_A = sum(Z_mat(:,5))+F_transf;
    Z_W = sum(Z_mat(:,6))+F_W;
    Z_T = sum(Z_mat(:,7))+F_T;
    disp('cost items: ')
    [Z_L,Z_V,Z_M,Z_A,Z_W,Z_T]
    disp('total cost: ')
    Z = Z_L + Z_V + Z_M + Z_A + Z_W + Z_T

    scen = f;
    run P_plots.m
    
    run M2_Outputs.m

    filename = strcat(num2str(rho_0),'-',num2str(gamma),'-',num2str(R),'-',num2str(v_p));
    switch scen
        case 0
            save(strcat(filename,'-','onlyBackbone.mat'))
        case 1
            save(strcat(filename,'-','fixedFMLM.mat'))
        case 2
            save(strcat(filename,'-','demrespFMLM.mat'))
        case 3
            save(strcat(filename,'-','adaptive.mat'))
        otherwise
            error('invalid "f" value');
    end
    
end %end "for" cycle 

toc % end timer