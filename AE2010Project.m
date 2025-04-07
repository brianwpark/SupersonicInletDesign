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
M_inlet_exit_suprsonic = 0.55;
M_vehicle_supersonic = 1.60;

%Baseline Design Parameters
Turn_angle_actual = deg2rad(5.5); % degrees




%BASELINE CALCULATION
%lift off calculation
%Area calculation
%speed of sound calculation
a_entrance_liftoff = sqrt(gamma*R_Specific*T_sea); % m/s
%Velocity calculation
V_liftoff = M_inlet_entrance_liftoff*a_entrance_liftoff; % m/s
%density calculation
rho = P_sea/(R_Specific*T_sea); % kg/m^3
%Area calculation
A_entrance = m_dot_liftoff/(rho*V_liftoff); % m^2
%Height and Width calculation
S = sqrt(A_entrance/2.5); % m
H1 = 2.5*S; % m

%Horiznotal Length calculation
Turn_angle_rad = deg2rad(5.5); % radians
L = H1*tan(Turn_angle_rad); % m


%Exit height calculation (currently thinking)
T0_exit_liftoff = T_sea*(1+(gamma-1)/2*M_inlet_entrance_liftoff^2); % K
a_exit_liftoff = sqrt(gamma*R_Specific*T0_exit_liftoff);
V_exit_liftoff = M_vehicle_liftoff*a_exit_liftoff; % m/s
p0_exit_liftoff = P_sea*(1+(gamma-1)/2*M_inlet_entrance_liftoff^2)^(gamma/(gamma-1)); % Pa
rho_exit_liftoff = p0_exit_liftoff/(R_Specific*T0_exit_liftoff); % kg/m^3
A_exit = m_dot_liftoff/(rho_exit_liftoff*V_exit_liftoff); % m
H2 = A_exit/S; % m

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





% BASELINE OFF DESIGN CALCULATION
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

%Density calculation
rho_ambient = P_ambient/(R_Specific*T_ambient); % kg/m^3

%oblique shock relations
for i = 1:length(Shock_formation)
    M_1_n = M_1(i)*sin(theta);
    M_1_a = M_1(i)
    if Shock_formation(i) == 0 % Detached Shock
        % 2nd region basic calculation
        [T_2,P_2,M_2] = NormalShockCalc(T_ambient,P_ambient,rho_ambient,M_1_n,gamma);
    elseif Shock_formation(i) == 1 % Attached Shock
        % 2nd region basic calculation
        [T_2,P_2,rho_2,M_2_n] = NormalShockCalc(T_ambient,P_ambient,rho_ambient,M_1_n,gamma);
        M_2 = M_2_n/sin(theta-Turn_angle_actual);
    end

    %3rd region basic calculation
    [T_3,P_3,rho_3,M_3] = NormalShockCalc(T_2,P_2,rho_2,M_2,gamma);
        
    %Exit nozzle condition calculation
    T0_3(i) = T_3*(1+(gamma-1)/2*M_3^2);
    P0_3(i) = P_3*(1+(gamma-1)/2*M_3^2)^(gamma/(gamma-1));

    V_entrance = sqrt(R_Specific*T_2*gamma)*M_2; % m/s
    m_dot_inlet(i) = rho_3*V_entrance*A_entrance;

    P0_ambient = P_ambient*(1+(gamma-1)/2*M_vehicle_supersonic^2)^(gamma/(gamma-1)); % Pa
    P_loss_fraction(i) = P0_3(i)/P0_ambient; % Pa
end
% Combined plot of T0_3 and P0_3
figure;
yyaxis left
plot(M_1, T0_3, 'b-', 'LineWidth', 2);
ylabel('Stagnation Temperature T0_3 (K)');

yyaxis right
plot(M_1, P0_3, 'r-', 'LineWidth', 2);
ylabel('Stagnation Pressure P0_3 (Pa)');

grid on;
xlabel('Mach Number (M_1)');
title('Stagnation Temperature and Pressure vs. Mach Number');
legend('Temperature', 'Pressure', 'Location', 'best');

% Plot M_1 vs Fractional Stagnation Pressure Loss
figure;
plot(M_1, P_loss_fraction, 'g-', 'LineWidth', 2);
grid on;
xlabel('Mach Number (M_1)');
ylabel('Fractional Stagnation Pressure Loss (P0_3/P0_{ambient})');
title('Fractional Stagnation Pressure Loss vs. Mach Number');

% Plot M_1 vs Mass Flow Rate
figure;
plot(M_1, m_dot_inlet, 'm-', 'LineWidth', 2);
grid on;
xlabel('Mach Number (M_1)');
ylabel('Mass Flow Rate (kg/s)');
title('Inlet Mass Flow Rate vs. Mach Number');

function [T_2,P_2,rho_2,M_2] = NormalShockCalc(T_1,P_1,rho_1,M_1,gamma)
    % This function calculates the properties of a normal shock wave.
    % Calculating M2, T2, P2, and rho2 from M1, T1, and P1.
    M_2 = sqrt((M_1^2+2/(gamma-1))/((2*gamma)/(gamma-1)*M_1^2-1));
    T_2 = T_1*(1+((gamma-1)/2)*M_1^2)/(1+((gamma-1)/2)*M_2^2);
    P_2 = P_1*(1+gamma*M_1^2)/(1+gamma*M_2^2);
    rho_2 = (M_1/M_2)*sqrt(T_1/T_2);
end