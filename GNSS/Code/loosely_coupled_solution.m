clc
clear
close

% add function and data that i  s given by paul.
addpath("Functions_Given", "Data_Given");
Define_Constants;
% Dataload
Pseudo_ranges_data = readmatrix("Pseudo_ranges.csv");
Pseudo_ranges_rates_data = readmatrix("Pseudo_range_rates.csv");
DR_Data = readmatrix("Dead_reckoning.csv");

% calculate the GNSS Solution, 
% return t L_b, lambda_b, h_b, v_N, v_E, Psi
GNSS_Solution = GNSS_Solver(Pseudo_ranges_data, Pseudo_ranges_rates_data);

% calculate the heading Solution, use gyro and meg ,return psi
psi_C = Gyro_Mag_KF_Solver(DR_Data);

% calculate the Dead reconing Solution, 
% return t L_b, lambda_b, h_b, v_N, v_E, Psi
DR_Solution = DR_Solver(DR_Data, GNSS_Solution, psi_C);

% couple the two solutions loosely
LCKF_Solution = LC_KF_Inte_Solver(DR_Solution, GNSS_Solution);

% couple the two solutions loosely with close loop
LCKF_Close_Solution = LC_KF_Inte_Close_Solver(DR_Data, GNSS_Solution, psi_C);

% plot figure
color = ['r', 'k', 'b', 'g'];
Result_Plot('GNSS', GNSS_Solution.L_b*rad_to_deg, ...
             GNSS_Solution.lambda_b*rad_to_deg, color(1));
Result_Plot('DR', DR_Solution.L_b*rad_to_deg, ...
             DR_Solution.lambda_b*rad_to_deg, color(2));
Result_Plot('Loosely Coupled Kalman Filter', ...
             LCKF_Solution.L_b*rad_to_deg, ...
             LCKF_Solution.lambda_b*rad_to_deg, color(3));
Result_Plot('Loosely Coupled closed loop Kalman Filter', ...
             LCKF_Close_Solution.L_b*rad_to_deg, ...
             LCKF_Close_Solution.lambda_b*rad_to_deg, color(4));
legend("DR with correction",'GNSS', 'DR', 'Loosely Coupled Kalman Filter', ...
       "L-C closed-loop Kalman Filter")

% store file
opres = [LCKF_Solution.t,LCKF_Solution.L_b'*rad_to_deg, ...
         LCKF_Solution.lambda_b'*rad_to_deg,LCKF_Solution.v_N', ...
         LCKF_Solution.v_E',psi_C'*rad_to_deg];
writematrix(opres,"Openloop_Output_Profile.csv");
cpres = [LCKF_Close_Solution.t,LCKF_Close_Solution.L_b'*rad_to_deg, ...
         LCKF_Close_Solution.lambda_b'*rad_to_deg,LCKF_Close_Solution.v_N', ...
         LCKF_Close_Solution.v_E',psi_C'*rad_to_deg];
writematrix(cpres,"Closedloop_Output_Profile.csv");



