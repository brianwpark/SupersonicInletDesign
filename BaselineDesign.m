%Design Requirements
gamma = 1.385;
M_mol = 29.1;
M_mol_KG = M_mol/1000; % kg/mol
R = 8.314; % J/(mol*K)
R_Specific = R/M_mol_KG; % J/(kg*K)

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
M_vehicle_supersonic = 1.60;

%Baseline Design Parameters
Turn_angle = 5.5; % degrees

%lift off calculation
%Area calculation
%speed of sound calculation
a_liftoff = sqrt(gamma*R_Specific*T_sea); % m/s
%Velocity calculation
V_liftoff = M_inlet_entrance_liftoff*a_liftoff; % m/s
%density calculation
rho = P_sea/(R_Specific*T_sea); % kg/m^3
%Area calculation
A_entrance = m_dot/(rho*V_liftoff); % m^2
%Height and Width calculation
S = sqrt(A_entrance/2.5); % m
H1 = 2.5*S; % m

%Horiznotal Length calculation
Turn_angle_rad = deg2rad(5.5); % radians
L = H1*tan(Turn_angle_rad); % m

%Exit height calculation (currently thinking)
H2 = 1.0; %dummy value, to be calculated later

%Supersonic inlet calculation
%Mass flow rate calculation
a_supersonic = sqrt(gamma*R_Specific*T_ambient); % m/s
V_Supersonic = M_inlet_exit_suprsonic*a_supersonic; % m/s
Rho_supersonic = P_ambient/(R_Specific*T_ambient); % kg/m^3
m_dot_supersonic = Rho_supersonic*V_Supersonic*A_entrance; % kg/s
%Stagnation pressure calculation
P2_supersonic = (((2*gamma)/(gamma+1))*(M_vehicle_supersonic^2)*(sin(Turn_angle_rad)^2)-((gamma-1)/(gamma+1)))*P_ambient; % Pa
P0_supersonic = (1+(gamma-1)/2*(M_vehicle_supersonic^2))^(gamma/(gamma-1))*P_ambient; % Pa

% Create a table with the key parameters
paramNames = {'Inlet Height (H1)', 'Width (S)', 'Length (L)', 'Exit Height (H2)', ...
              'Stagnation Pressure (P0)', 'Mass Flow Rate'};
paramUnits = {'m', 'm', 'm', 'm', 'Pa', 'kg/s'};
paramValues = [H1, S, L, H2, P0_supersonic, m_dot_supersonic];

designTable = table(paramNames', paramUnits', paramValues', ...
                   'VariableNames', {'Parameter', 'Unit', 'Value'});
disp(designTable);