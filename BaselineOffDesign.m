%{
    This code is for calculating stagnation temperature and pressure, fractional stagnation pressure loss, and mass flow rate entering the inlet.
%}



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
Turn_angle_actual = deg2rad(5.5); % rad

%Detached or Attached Shock
%https://math.stackexchange.com/questions/4699275/finding-maximum-deflection-angle-for-oblique-shock-waves
syms shock_wave_angle_rad
Turn_angle = atan(2*cot(shock_wave_angle_rad)*(((M_vehicle_supersonic)^2*sin(shock_wave_angle_rad)^2-1)/((M_vehicle_supersonic)^2*(gamma+cos(2*shock_wave_angle_rad))))); %rad
d_turn_d_shock = diff(Turn_angle, shock_wave_angle_rad);
d_turn_d_shock_func = matlabFunction(d_turn_d_shock); % Convert symbolic function to numeric function handle
Solved_shock_angle = fzero(d_turn_d_shock_func, [0.1, pi/2]);
Max_turn_angle = atan(2*cot(Solved_shock_angle)*(((M_vehicle_supersonic)^2*sin(Solved_shock_angle)^2-1)/((M_vehicle_supersonic)^2*(gamma+cos(2*Solved_shock_angle)))));

if Turn_angle_actual > Max_turn_angle
    disp('The oblique shock is a Detached Shock')
else
    disp('The oblique shock is a Attached Shock')
end