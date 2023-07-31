function [x_pos,latitude,longitude,height,velocity,H_G,delta_z] = ...
                GNSS_LS_Solver(sat_names,pseudo_ranges_epoch, ...
                               pseudo_range_rates_epoch,x_pri)
    % Constant
    r_eb_e = x_pri(1:3);    
    v_eb_e = x_pri(4:6);
    delta_rho = x_pri(7);
    dot_delta_rho = x_pri(8);
    time = pseudo_ranges_epoch(1);
    num_sat = length(pseudo_ranges_epoch)-1;
    omega_ie = 7.292115e-5;
    Omega_ie = Skew_symmetric([0,0,omega_ie]);
    c = 299792458;
    max_it = 10;
    R_a = zeros(1,8);
    for i = 1:max_it
        % b)
        Sat_r_es_e = zeros(3,num_sat);
        Sat_v_es_e = zeros(3,num_sat);   
        for j = 2:num_sat+1   
            [sat_r_es_e,sat_v_es_e] = Satellite_position_and_velocity(time,sat_names(j));
            Sat_r_es_e(:,j-1) = sat_r_es_e';
            Sat_v_es_e(:,j-1) = sat_v_es_e';
        end

        % c) & d)
        C_Ie = eye(3);
        U_a = zeros(3,num_sat);
        dot_R_a = zeros(1,num_sat);
        for j = 1:length(Sat_r_es_e)
            if i == 1
                R_a(j) = sqrt((C_Ie*Sat_r_es_e(:,j)-r_eb_e)'* ...
                             (C_Ie*Sat_r_es_e(:,j)-r_eb_e));
                C_Ie = [1, omega_ie*R_a(j)/c, 0;
                    -omega_ie*R_a(j)/c, 1, 0;
                    0, 0, 1];
            else
                C_Ie = [1, omega_ie*R_a(j)/c, 0;
                    -omega_ie*R_a(j)/c, 1, 0;
                    0, 0, 1];
            end
            R_a(j) = sqrt((C_Ie*Sat_r_es_e(:,j)-r_eb_e)'* ...
                         (C_Ie*Sat_r_es_e(:,j)-r_eb_e));
            U_a(:,j) = (C_Ie*Sat_r_es_e(:,j)-r_eb_e)/R_a(j);
            dot_R_a(j) = U_a(:,j)'*(C_Ie*(Sat_v_es_e(:,j)+ ...
                         Omega_ie*Sat_r_es_e(:,j))- ...
                         (v_eb_e+Omega_ie*r_eb_e));
        end
        
        
        % e)
        x_hat = [r_eb_e; delta_rho];
        dot_x_hat = [v_eb_e; dot_delta_rho];
        delta_z = zeros(num_sat,1);
        dot_delta_z = zeros(num_sat,1);
        H_G = zeros(num_sat,4);
        for j = 2:num_sat+1
            delta_z(j-1) =  pseudo_ranges_epoch(j)-R_a(j-1)-delta_rho;
            dot_delta_z(j-1) =  pseudo_range_rates_epoch(j)-dot_R_a(j-1)-dot_delta_rho;
            H_G(j-1,:) = [-U_a(:,j-1)',1]; 
        end
        
        % f)
        x_hat = x_hat + (H_G'*H_G)\H_G'*delta_z;
        dot_x_hat = dot_x_hat + (H_G'*H_G)\H_G'*dot_delta_z;
        
        % g)
        [latitude,longitude,height,velocity] = pv_ECEF_to_NED(x_hat(1:3),dot_x_hat(1:3));
        
        % update
        x_pos = [x_hat(1:3);dot_x_hat(1:3);x_hat(4);dot_x_hat(4)];
        if norm(r_eb_e-x_hat(1:3)) < 0.1
            break;
        end
        r_eb_e = x_hat(1:3);
        v_eb_e = dot_x_hat(1:3);
        delta_rho = x_hat(4);
        dot_delta_rho = dot_x_hat(4);
    end
end



