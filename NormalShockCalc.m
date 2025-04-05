function [T_2,rho_2,P_2,M_2] = NormalShockCalc(T_1,rho_1,P_1,M_1,gamma)
    % This function calculates the properties of a normal shock wave.
    % Calculating M2, T2, P2, and rho2 from M1, T1, P1, and rho1.
    M_2 = (M_1^2+2/(gamma-1))/((2*gamma)/(gamma-1)*M_1^2-1);
    T_2 = T_1*(1+((gamma-1)/2)*M_1^2)/(gamma*M_2^2-(gamma-1)/2);
    P_2 = (1+gamma*M_1^2)/(1+gamma*M_2^2);
    rho_2 = (M_1/M_2)*sqrt(T_1/T_2)*rho_1;
end