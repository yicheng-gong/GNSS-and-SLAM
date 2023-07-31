function LCKF_Solution = LC_KF_Inte_Close_Solver(DR_Data, GNSS_Solution, Psi_C)
    Define_Constants;
    % init DR data param to loop
    h_0 = GNSS_Solution.h_b(1); % check here for every GNSS data
    v_bar(1) = (DR_Data(1,4) + DR_Data(1,5))/2 ;
    t = DR_Data(:,1);
    L_k_D(1) = GNSS_Solution.L_b(1);
    lambda_k_D(1) = GNSS_Solution.lambda_b(1);
    v_N_D(1) = v_bar(1) * cos(psi(1));
    v_E_D(1) = v_bar(1) * sin(psi(1));
    % load GNSS data 
    L_k_G = GNSS_Solution.L_b;
    lambda_k_G = GNSS_Solution.lambda_b;
    h_k_G = GNSS_Solution.h_b;
    v_N_G= GNSS_Solution.v_N;
    v_E_G= GNSS_Solution.v_E;

    % init state and covariant 
    [R_N, R_E]= Radii_of_curvature(L_k_G(1));
    x_pos = zeros(4,1); % v_N, v_E, L, lambda
    P_pos = [eye(2)*0.1^2, zeros(2);
             zeros(2), [(10/(R_N+h_k_G(1)))^2, 0;
                        0, (10/((R_E+h_k_G(1)) * cos(L_k_G(1))))^2]];
    tau_s = 0.5;
    S_DR = 0.2;
    sigma_Gr = 5;
    sigma_Gv = 0.02;

    %%%%%%% prepare kalman filter parameters
    Phi = eye(4);
    Phi(3,1) = 0.5/(R_N + h_k_G(1));
    Phi(4,2) = 0.5/((R_E + h_k_G(1)) * cos(L_k_D(1)));
    % define system noise covariance % 
    Q = [S_DR*tau_s, 0 , 0.5*S_DR*tau_s^2/(R_N+h_k_G(1)), 0;
         0, S_DR*tau_s, 0, 0.5*S_DR*tau_s^2/((R_E+h_k_G(1))*cos(L_k_D(1)));
         0.5*S_DR*tau_s^2/(R_N+h_k_G(1)), 0, 1/3*S_DR*tau_s^3/(R_N+h_k_G(1))^2, 0;
         0, 0.5*S_DR*tau_s^2/((R_E+h_k_G(1))*cos(L_k_D(1))), 0, ...
         1/3*S_DR*tau_s^3/((R_E+h_k_G(1))^2*cos(L_k_D(1))^2)];
    
    % Propagate the state estimates and the error covariance matrix
    x_pri = Phi*x_pos;
    P_pri = Phi*P_pos*Phi' + Q;

    % Compute the measurement matrix
    H = [0 0 -1 0;
         0 0 0 -1;
         -1 0 0 0;
         0 -1 0 0];
    
    % Compute the measurement noise covariance
    R = [sigma_Gr^2/(R_N+ h_k_G(1))^2, 0, 0, 0;
        0, sigma_Gr^2 / (R_E+h_k_G(1))^2 / cos(L_k_G(1))^2, 0, 0;
        0, 0,  sigma_Gv^2 ,0;
        0, 0, 0, sigma_Gv^2];
    K = P_pri*H'*inv(H*P_pri*H' + R);
        
    % Formulate the measurement innovation
    delta_z = [L_k_G(1) - L_k_D(1);
               lambda_k_G(1) - lambda_k_D(1);
               v_N_G(1) - v_N_D(1);
               v_E_G(1) - v_E_D(1);
               ] - H*x_pri;

    % Update the state estimates, the error covariance matrix
    x_pos = x_pri + K*delta_z;
    P_pos = (eye(4) - K*H)*P_pri;
    v_N_C(1) = v_N_D(1) - x_pos(1);
    v_E_C(1) = v_E_D(1) - x_pos(2);
    L_k_C(1) = L_k_D(1)- x_pos(3);
    lambda_k_C (1) = lambda_k_D(1) - x_pos(4);

    %%%%%%% start kalman filter 
    for k = 2:length(t)
        % calculate DR solution for this epoch
        v_bar(k) = (DR_Data(k,4) + DR_Data(k,5))/2 ;
        v_N_bar = 0.5*(cos(Psi_C(k)) + cos(Psi_C(k-1)))*v_bar(k);
        v_E_bar = 0.5*(sin(Psi_C(k)) + sin(Psi_C(k-1)))*v_bar(k);
        [R_N,R_E] = Radii_of_curvature(L_k_C(k-1));
        L_k_D(k) = L_k_D(k-1) + (v_N_bar*(t(k) - t(k-1))) / (R_N + h_0);
        lambda_k_D(k) = lambda_k_D(k-1) + (v_E_bar*(t(k) - t(k-1)) / ...
                        ((R_E + h_0) * cos(L_k_D(k)))); 
        v_N_D(k) = 1.7 * v_N_bar - 0.7 * v_N_D(k-1);   
        v_E_D(k) = 1.7 * v_E_bar - 0.7 * v_E_D(k-1);

        % calculate transition matrix 
        [R_N, R_E]= Radii_of_curvature(L_k_G(k));
        Phi = eye(4);
        Phi(3,1) = 0.5/(R_N + h_k_G(k-1));
        Phi(4,2) = 0.5/((R_E + h_k_G(k-1)) * cos(L_k_D(k-1)));

        % define system noise covariance 
        Q = [S_DR*tau_s, 0 , 0.5*S_DR*tau_s^2/(R_N+h_k_G(k-1)), 0;
             0, S_DR*tau_s, 0, 0.5*S_DR*tau_s^2/((R_E+h_k_G(k-1))*cos(L_k_D(k-1)));
             0.5*S_DR*tau_s^2/(R_N+h_k_G(k-1)), 0, 1/3*S_DR*tau_s^3/(R_N+h_k_G(k-1))^2, 0;
             0, 0.5*S_DR*tau_s^2/((R_E+h_k_G(k-1))*cos(L_k_D(k-1))), 0, ...
             1/3*S_DR*tau_s^3/((R_E+h_k_G(k-1))^2*cos(L_k_D(k-1))^2)];
        
        % Propagate the state estimates and the error covariance matrix
        x_pri = Phi*x_pos;
        P_pri = Phi*P_pos*Phi' + Q;
        
        % Compute the measurement matrix
        H = [0 0 -1 0;
             0 0 0 -1;
             -1 0 0 0;
             0 -1 0 0];
        % Compute the measurement noise covariance
        R = [sigma_Gr^2/(R_N+ h_k_G(k))^2, 0, 0, 0;
            0, sigma_Gr^2 / (R_E+h_k_G(k))^2 / cos(L_k_G(k))^2, 0, 0;
            0, 0,  sigma_Gv^2 ,0;
            0, 0, 0, sigma_Gv^2];
        
        % Compute the Kalman gain matrix
        K = P_pri*H'*inv(H*P_pri*H' + R);
        
        % Formulate the measurement innovation
        delta_z = [L_k_G(k) - L_k_D(k);
                   lambda_k_G(k) - lambda_k_D(k);
                   v_N_G(k) - v_N_D(k);
                   v_E_G(k) - v_E_D(k);
                   ] - H*x_pri;

        % Update the state estimates, the error covariance matrix
        x_pos = x_pri + K*delta_z;
        P_pos = (eye(4) - K*H)*P_pri;
        v_N_C(k) = v_N_D(k) - x_pos(1);
        v_E_C(k) = v_E_D(k) - x_pos(2);
        L_k_C(k) = L_k_D(k)- x_pos(3);
        lambda_k_C(k) = lambda_k_D(k) - x_pos(4);
        
        if mod(k,5) == 0
            x_pos = [0;0;0;0];
            v_N_D(k) = v_N_C(k);
            v_E_D(k) = v_E_C(k);
            L_k_D(k) = L_k_C(k);
            lambda_k_D(k) = lambda_k_C(k);
        end
    end
    figure
    plot(lambda_k_D*rad_to_deg, L_k_D*rad_to_deg, '-');
    title("DR Result after Closed-loop correction", ' Position');
    xlabel('Latitude (deg)');
    ylabel('Longitude (deg)');
    grid on;
    hold on
    LCKF_Solution.v_N = v_N_C;
    LCKF_Solution.v_E = v_E_C;
    LCKF_Solution.L_b = L_k_C;
    LCKF_Solution.lambda_b = lambda_k_C;
    LCKF_Solution.t = t;
end