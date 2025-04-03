%Design Requirements
gamma = 1.385;
M_mol = 29.1;
M_mol_KG = Mmol/1000; % kg/mol
R = 8.314; % J/(mol*K)
R_Specific = R/MmolKG; % J/(kg*K)

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

%Detached or Attached Shock
%https://math.stackexchange.com/questions/4699275/finding-maximum-deflection-angle-for-oblique-shock-waves

if Turn_angle > Max_turn_angle
    disp('The oblique shock is a Detached Shock')
else
    disp('The oblique shock is a Attached Shock')
end