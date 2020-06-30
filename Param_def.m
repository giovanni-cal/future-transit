%% Parameter Definition
%% Conversions
km_mi = 1.61; % km/miles

%% Constants
%[ rail metro bus drt ]
v_p   =  3.6         ; % [km/h]
v     = [70 70 30 30]; % [km/h] T = d/v +n_s*tau_s+n_pax*tau_pax
tau_s = [0.025 0.0167 0.0083 0.0056] ; % [h] (90 sec - rail; 60 sec - metro; 30 sec bus; 20 sec drt)
tau_p =  0.00028     ; % [h] ( 1 sec)
d     =  0.3           ; % [km]

%cost coefficients - inspired by Chen and Daganzo (2015)

c_L =  [ 200 600 10 10 0]; %[eur/km-h] %infrastructure
c_V =  [ 6 6 1 1 0]; %[eur/veh-km]]
c_M =  [150 150 30 30 0]; %[eur/veh-h]
Cap =  [ 1500 1500 100 30]; %[pax/veh]

% Value of Time scenarios (5)
mu_T = 10; %[eur/h] ... [ 5:5:25 ] in-vehicle
mu_W = 1.5 * mu_T; %[eur/h] waiting
mu_A = 1.5 * mu_T  ; %[eur/h] access

% area scenarios (10)
R = 25; %[km] ... [ 25:5:25 ]

% demand density scenarios (20)
rho_0 = 300; %[trips/km2-h] ... [ 30:90:300 ]
gamma = 0.12;
% rho_0 = [ 50:50:1000 ;
%           25:25:500 ;
%           10:10:200 ;
%           5:5:100 ] ;
%%%demand = 2 * pi * rho_0 * (1 - exp(-gamma*R)*(gamma*R + 1)) / gamma^2
Dem_x = @(x) (2*pi*x*rho_0.*exp(-gamma*x));
tot_Dem = integral(Dem_x,0,R)



