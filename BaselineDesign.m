%Design Requirements
gamma = 1.385;
M_mol = 29.1;
M_mol_KG = Mmol/1000; % kg/mol
R = 8.314; % J/(mol*K)
R_Specific = R/MmolKG; % J/(kg*K)

Turn_angle = 5.5; % degrees

%Liftoff conditions
T_sea = 288.15; % K
P_sea = 101325; % Pa
M_inlet_entrance_liftoff = 0.55;
M_vehicle_liftoff = 0.31;
m_dot = 63.6; % kg/s

%Supersonic cruise conditions
P_ambient = 11.6e3; % Pa
T_ambient = 216e3; % K
M_inlet_exit_suprsonic = 0.55;

%Baseline Design Parameters
Turn_angle = 5.5; % degrees

%lift off calculation
%Area calculation
%speed of sound calculation
a = sqrt(gamma*R_Specific*T_sea); % m/s
%Velocity calculation
v = M_inlet_liftoff*a; % m/s
%density calculation
rho = P_sea/(R_Specific*T_sea); % kg/m^3
%Area calculation
A = m_dot/(rho*v); % m^2
%Height and Width calculation
S = sqrt(A/2.5); % m
H1 = 2.5*S; % m

%Horiznotal Length calculation
Turn_angle_rad = deg2rad(5.5); % radians
L = H1*tan(Turn_angle_rad); % m

%Exit height calculation (currently thinking)

%Supersonic inlet calculation
%Mass flow rate calculation
a = sqrt(gamma*R_Specific*T_ambient); % m/s