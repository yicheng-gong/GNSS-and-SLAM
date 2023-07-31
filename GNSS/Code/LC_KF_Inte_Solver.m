function LCKF_Solution = LC_KF_Inte_Solver(DR_Solution, GNSS_Solution)
    
    L_k_D = DR_Solution.L_b;
    lambda_k_D = DR_Solution.lambda_b;
    v_N_D = DR_Solution.v_N;
    v_E_D = DR_Solution.v_E;
    
    t = DR_Solution.t;
    L_k_G = GNSS_Solution.L_b;
    lambda_k_G = GNSS_Solution.lambda_b;
    h_k_G = GNSS_Solution.h_b;
    v_N_G= GNSS_Solution.v_N;
    v_E_G= GNSS_Solution.v_E;

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
    [R_N, R_E]= Radii_of_curvature(L_k_G(1));
    Phi = eye(4);
    Phi(3,1) = 0.5/(R_N + h_k_G(1));
    Phi(4,2) = 0.5/((R_E + h_k_G(1)) * cos(L_k_D(1)));
    % define system noise covariance % 
    Q = [S_DR*tau_s, 0 , 0.5*S_DR*tau_s^2/(R_N+h_k_G(1)), 0;
         0, S_DR*tau_s, 0, 0.5*S_DR*tau_s^2/((R_E+h_k_G(1))*cos(L_k_D(1)));
         0.5*S_DR*tau_s^2/(R_N+h_k_G(1)), 0, 1/3*S_DR*tau_s^3/(R_N+h_k_G(1))^2, 0;
         0, 0.5*S_DR*tau_s^2/((R_E+h_k_G(1))*cos(L_k_D(1))), ...
         0, 1/3*S_DR*tau_s^3/((R_E+h_k_G(1))^2*cos(L_k_D(1))^2)];
    
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
        [R_N, R_E]= Radii_of_curvature(L_k_G(k));
        Phi = eye(4);
        Phi(3,1) = 0.5/(R_N + h_k_G(k-1));
        Phi(4,2) = 0.5/((R_E + h_k_G(k-1)) * cos(L_k_D(k-1)));
        % define system noise covariance % 
        Q = [S_DR*tau_s, 0 , 0.5*S_DR*tau_s^2/(R_N+h_k_G(k-1)), 0;
             0, S_DR*tau_s, 0, 0.5*S_DR*tau_s^2/((R_E+h_k_G(k-1))*cos(L_k_D(k-1)));
             0.5*S_DR*tau_s^2/(R_N+h_k_G(k-1)), 0, 1/3*S_DR*tau_s^3/(R_N+h_k_G(k-1))^2, 0;
             0, 0.5*S_DR*tau_s^2/((R_E+h_k_G(k-1))*cos(L_k_D(k-1))), ...
             0, 1/3*S_DR*tau_s^3/((R_E+h_k_G(k-1))^2*cos(L_k_D(k-1))^2)];
        
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
        lambda_k_C (k) = lambda_k_D(k) - x_pos(4);
        
    end
    LCKF_Solution.v_N = v_N_C;
    LCKF_Solution.v_E = v_E_C;
    LCKF_Solution.L_b = L_k_C;
    LCKF_Solution.lambda_b = lambda_k_C;
    LCKF_Solution.t = t;
end