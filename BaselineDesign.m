%Design Requirements
gamma = 1.385;
M_mol = 29.1;
M_mol_KG = Mmol/1000; % kg/mol
R = 8.314; % J/(mol*K)
R_Specific = R/MmolKG; % J/(kg*K)
M_inlet_liftoff = 0.55;

Turn_angle = 5.5; % degrees

%Lift off parameters
T_sea = 288.15; % K
P_sea = 101325; % Pa

%required parameters
m_dot = 63.6; % kg/s

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


%Horiznotal Length calculation
Turn_angle_rad = deg2rad(5.5); % radians
L = H1*tan(Turn_angle_rad); % m
%Exit height calculation

