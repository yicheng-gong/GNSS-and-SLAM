function psi_C = Gyro_Mag_KF_Solver(DR_Data)
    Define_Constants;
    omega_G = DR_Data(:,6);
    psi_M = DR_Data(:,7)*deg_to_rad;
    psi_G = zeros(length(psi_M), 1);
    psi_G(1) = psi_M(1);
    Tau_s = 0.5;
    for i = 2:length(psi_M)
        psi_G(i) = psi_G(i-1) + omega_G(i-1) * Tau_s;
    end
    
    % assume state x = [delta_psi, b]  b = 1deg/s
    x_pos = [0; deg_to_rad];
    P_pos = [0.01^2,  0;
             0, 0.001^2];
    % start kalman 
    for k = 1:length(psi_G)
        S_rg = 1e-4;
        S_bgd = deg_to_rad;
        Phi = [1, Tau_s;
               0, 1];
        Q = [S_rg*Tau_s + 1/3*S_bgd*Tau_s^3, 0.5*S_bgd*Tau_s^2;
            0.5*S_bgd*Tau_s^2, S_bgd *Tau_s];

        % transite the model
        x_pri = Phi*x_pos;
        P_pri = Phi*P_pos*Phi' + Q;

        % define H and R
        H = [-1 0];
        R = (4 * deg_to_rad)^2;
        K = P_pri*H'/(H*P_pri*H' + R);

        % calculate measurement inovation
        delta_z = (psi_M(k) - psi_G(k)) - H * x_pri;

        % calculate posterir
        x_pos = x_pri + K*delta_z;
        P_pos = (eye(2) - K*H)*P_pri;
        
        % get corrected value
        psi_C(k) = psi_G(k) - x_pos(1);
    end

    figure
    plot(DR_Data(:,1),psi_M, Color='r');
    hold on
    plot(DR_Data(:,1),psi_G, Color='k');
    plot(DR_Data(:,1),psi_C, Color='b');
    xlabel('time (sec)')
    ylabel('Heading (rad)')
    title('Heading Calculation')
    legend("meg","gyro","correct", Location="southwest")
end