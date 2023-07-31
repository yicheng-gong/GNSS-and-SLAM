function GNSS_Solution = GNSS_Solver(Pseudo_ranges_data, Pseudo_ranges_rates_data)
    % Pseudo_ranges_data:
    % col 2-9, row 2-end : satellete data
    % col 1, row 2-end : time

    % Pseudo_ranges_rate_data:
    % col 2-9, row 2-end : satellete data
    % col 1, row 2-end : time
    Define_Constants;

    num_epoch = size(Pseudo_ranges_data,1);
    num_epoch = num_epoch -1;
    % Task2A (a) Initialise the Kalman filter
    [x_pos,P_pos] = Initialise_GNSS_KF(Pseudo_ranges_data, ...
                                       Pseudo_ranges_rates_data);
    
    for epoch =  2:num_epoch+1
        % Outlier detect and remove
        [satellite_name_list, pseudo_ranges, pseudo_range_rates ] = ...
                   Detect_Outlier_And_Remove(x_pos,epoch, ...
                             Pseudo_ranges_data, Pseudo_ranges_rates_data);
        num_satellite = length(satellite_name_list);
        % Task2A (b) Compute the transition matrix
        tau_s = 0.5;
        Phi = [eye(3), tau_s*eye(3), zeros(3,2);
               zeros(3), eye(3), zeros(3,2);
               zeros(1,6), 1, tau_s;
               zeros(1,7), 1];
        
        % Task2A (c) Compute the system noise covariance matrix
        S_a = 1;
        S_c_phi = 0.01;
        S_c_f = 0.04;
        Q = [1/3*S_a*tau_s^3*eye(3), 0.5*S_a*tau_s^2*eye(3), zeros(3,2);
             0.5*S_a*tau_s^2*eye(3), S_a*tau_s*eye(3), zeros(3,2);
             zeros(1,6), S_c_phi*tau_s + 1/3*S_c_f*tau_s^3, 0.5*S_c_f*tau_s^2;
             zeros(1,6), 0.5*S_c_f*tau_s^2, S_c_f*tau_s];
        
        % Task2A (d) Use the transition matrix to propagate the state estimates
        x_pri = Phi * x_pos;
        
        % Task2A (e) use this to propagate the error covariance matrix
        x_conv_pri = Phi*P_pos*Phi' + Q;
        
        % Task2A (f)(g) Predict the ranges from the approximate user position to each satellite
        % calculate satellite position
        r_ej = zeros(3,num_satellite);
        v_ej = zeros(3,num_satellite);
        for i =  1:num_satellite
            [position_satellite, velocity_satellite] = ...
                                 Satellite_position_and_velocity ...
                                     (Pseudo_ranges_data(epoch,1), ...
                                      satellite_name_list(i));
            r_ej(:,i) = position_satellite';
            v_ej(:,i) = velocity_satellite';
        end
        
        r_aj_pri = zeros(1,num_satellite);
        u_aj = zeros(3,num_satellite);
        r_aj_dot_pri = zeros(1,num_satellite);
        for i = 1:num_satellite
            % calculate C
%             if epoch == 2
                C = eye(3);
                r_aj_pri(i) = sqrt((C* r_ej(:,i) - x_pri(1:3))' * ...
                                   (C* r_ej(:,i) - x_pri(1:3)));
%             else
%                 C = eye(3);
%             end
            C(1,2) = omega_ie * r_aj_pri(i) /c;
            C(2,1) = -omega_ie * r_aj_pri(i) /c;
            r_aj_pri(i) = sqrt((C* r_ej(:,i) - x_pri(1:3))' * ...
                               (C* r_ej(:,i) - x_pri(1:3)));
            u_aj(:,i) = (C * r_ej(:,i) - x_pri(1:3)) / r_aj_pri(i);
            % Task2A (h) Predict the range rates
            r_aj_dot_pri(i) = u_aj(:,i)' * ...
                             (C * (v_ej(:,i) + Omega_ie*r_ej(:,i)) - ...
                              (x_pri(4:6) + Omega_ie* x_pri(1:3))) ;
        end
        
        % Task2A (i) Compute the measurement matrix
        H = [-u_aj', zeros(num_satellite,3), ones(num_satellite,1), zeros(num_satellite,1);
             zeros(num_satellite,3), -u_aj', zeros(num_satellite,1),ones(num_satellite,1)];
        
        % Task2A (j) Compute the measurement noise covariance
        sigma_rho = 10;
        sigma_r = 0.05;
        R_k = [eye(num_satellite)*sigma_rho^2 , zeros(num_satellite);
               zeros(num_satellite), eye(num_satellite)*sigma_r^2];
        
        % Task2A (k) Compute the Kalman gain
        K_k = x_conv_pri * H' / (H*x_conv_pri*H' + R_k);
        
        % Task2A (I) Formulate the measurement innovation
        delta_z = [(pseudo_ranges - r_aj_pri - x_pri(7))';
                   (pseudo_range_rates - r_aj_dot_pri - x_pri(8))'];
        
        % Task2A (m) Update the state estimates
        x_pos = x_pri + K_k * delta_z;
        
        % Task2A (n) Update the error covariance matrix
        P_pos = (eye(8) - K_k*H)*x_conv_pri;
        
        % Task2A (o) check solution
        [latitude(epoch-1), longitude(epoch-1) , ...
           height(epoch-1), v_ebn(:,epoch-1)] = pv_ECEF_to_NED(x_pos(1:3), x_pos(4:6));
        
    end
    GNSS_Solution.t = Pseudo_ranges_data(2:end,1);
    GNSS_Solution.L_b = latitude;
    GNSS_Solution.lambda_b = longitude;
    GNSS_Solution.h_b = height;
    GNSS_Solution.v_N = v_ebn(1,:);
    GNSS_Solution.v_E = v_ebn(2,:);

    figure
    plot(longitude*rad_to_deg, latitude*rad_to_deg, '-');
    title('GNSS Position Result');
    xlabel('Latitude (deg)');
    ylabel('Longitude (deg)');
    grid on;

    figure
    quiver(longitude*rad_to_deg,latitude*rad_to_deg, ...
                     GNSS_Solution.v_E,GNSS_Solution.v_N)
    grid on;
    xlabel('Longitude (deg)');
    ylabel('Latitude (deg)');
    title('GNSS Velocity Result')

end