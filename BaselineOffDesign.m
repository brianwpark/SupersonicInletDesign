%{
    This code is for calculating stagnation temperature and pressure, fractional stagnation pressure loss, and mass flow rate entering the inlet.
    This is how the Inlet is formed in supersonic flow:
    1st region, oblique shock wave, 2nd region, normal shock wave, 3rd region, exit nozzle.
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
m_dot_liftoff = 63.6; % kg/s

%Supersonic cruise conditions
P_ambient = 11.6e3; % Pa
T_ambient = 216; % K
% M_inlet_exit_suprsonic = 0.55;
M_vehicle_supersonic = 1.60;

%Baseline Design Parameters
Turn_angle_actual = deg2rad(5.5); % rad

%Detached or Attached Shock
%https://math.stackexchange.com/questions/4699275/finding-maximum-deflection-angle-for-oblique-shock-waves
M_1 = 1:0.1:1.8; % M_1 range for plotting
for i = 1:length(M_1)
    M_1_values = M_1(i);

    a = (M_1_values^2*(1+gamma)-4)/(2*M_1_values^2*gamma);
    b = (M_1_values^2*(gamma+1)+2)/(2*M_1_values^4*gamma);
    c = a+sqrt(a^2+4*b);


    Max_turn_angles(i) = atan(sqrt(2/c-1)*(M_1_values^2*c-2)/(M_1_values^2*(1+gamma-c)+2));
end

for i = 1:length(Max_turn_angles)
    Max_angle = Max_turn_angles(i);
    if Turn_angle_actual > Max_angle
        % fprintf('The oblique shock is a Detached Shock for %.2fth values\n',i);
        Shock_formation(i) = 1; % Detached Shock
    elseif Turn_angle_actual < Max_angle
        % fprintf('The oblique shock is a Attached Shock for %.2fth values\n',i);
        Shock_formation(i) = 0; % Attached Shock
    end
end

%Shock angle calculation
Lamda = sqrt((M_vehicle_supersonic^2-1)^2-3*(1+(gamma-1)/2*(M_vehicle_supersonic^2))*(1+(gamma+1)/2*(M_vehicle_supersonic^2))*(tan(Turn_angle_actual)^2));
x = (1/(Lamda^3))*((M_vehicle_supersonic^2-1)^3-9*(1+(gamma-1)/2*M_vehicle_supersonic^2)* (1+(gamma-1)/2*M_vehicle_supersonic^2+(gamma+1)/4*M_vehicle_supersonic^4)*(tan(Turn_angle_actual)^2));
alpha = 1;
theta = atan((M_vehicle_supersonic^2-1+2*Lamda*cos((4*pi*alpha+acos(x))/(3)))/(3*(1+(gamma-1)/(2)*M_vehicle_supersonic^2)*tan(Turn_angle_actual)));

%oblique shock relations
for i = 1:length(Shock_formation)
    M_1_n = M_1(i)*sin(theta);
    if Shock_formation(i) == 0 % Detached Shock
        % 2nd region basic calculation
        [T_2,P_2,M_2_n] = NormalShockCalc(T_ambient,P_ambient,M_1_n,gamma);
        M_2 = M_2_n/sin(theta-Turn_angle_actual);
        % Mass flow rate entering the inlet
        m_dot_inlet(i) = rho

        %3rd region basic calculation
        [T_3,P_3,M_3] = NormalShockCalc(T_2,P_2,M_2,gamma);
        
        %Exit nozzle condition calculation
        T0_3(i) = T_3*(1+(gamma-1)/2*M_3^2);
        P0_3(i) = P_3*(1+(gamma-1)/2*M_3^2)^(gamma/(gamma-1));

        m_dot_inlet(i) = rho_3
    elseif Shock_formation(i) == 1 % Attached Shock
        
    end
end