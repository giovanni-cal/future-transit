%% Constants Definition

%[ rail bus drt ]
v_p   = 4.5         ; % [km/h]
v     = [60 25 25]; % [km/h] T = D/v +n_s*tau_s+n_pax*tau_pax
tau_s = [0.0125 0.0083 0.0083] ; % [h] (45 sec rail & metro; 30 sec bus; 30 sec drt)
tau_p = 0.00056     ; % [h] ( 2 sec)
%d     = 0.3           ; % [km]
Wstr_max = 0.5; %km %mx width of a strip in the FMLM

%cost coefficients - inspired by Chen and Daganzo (2015)

c =  [ 600 10 10; 6 0.5 0.5; 120 50 50];
% infrastructure[eur/km-h]; operation[eur/veh-km]; fleet[eur/veh-h]
Cap =  [ 1200 80 80 ]; %[pax/veh]

% Value of Time scenarios (5)
mu = [15; 22.5; 30];
%in-vehicle[eur/h]; waiting[eur/h]; access[eur/h] 

% demand density scenarios (20)
%   rho_0 = (from input.txt); %[trips/km2-h]
%   gamma = (from input.txt);
%%% demand = 2 * pi * rho_0 * (1 - exp(-gamma*R)*(gamma*R + 1)) / gamma^2
Dem_x = @(x) (2*pi*x*rho_0.*exp(-gamma*x));
tot_Dem = integral(Dem_x,0,R)

%parameters for fmincon
%[Th_r,S_c,s,phi,H,h,d]
A = [-1,0,0,1,0,0,0;0,-1,1,0,0,0,0;0,0,-1,0,0,0,1]; b = [0;0;0];
%    phi <= Th_r ;  s <= S_c ;     d <= s;
A_eq = []; b_eq = []; %A_eq * DX = b_eq
lb = [0.02,0.5,0.5,0.02,0.042,0.033,0.2]; %lower bound of DX
ub = [pi,2.5,2.5,pi,0.5,0.5,2]; %upper bound of DX
options = optimoptions('fmincon','Display','off','MaxIterations',1000, ...
    'OptimalityTolerance',10^-6,'StepTolerance',10^-10,...
    'ConstraintTolerance',10^-6','UseParallel',1);
%'ConstraintTolerance',10^-6',OptimalityTolerance',10^-6,'StepTolerance',10^-10
%(Available algorithms: 'active-set', 'interior-point', 'sqp', 'sqp-legacy', 'trust-region-reflective')
options_glob = optimset('Display','off'); %for fminbnd



