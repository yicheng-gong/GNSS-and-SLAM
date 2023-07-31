function DR_Solution = DR_Solver(DR_Data, GNSS_Solution, psi_C)
    % DR_Data:
    % columns 1 time
    % columns 2-5: wheel-speed measurements
    % column 6: the gyroscope angular rate radians per second
    % column 7: contains the heading measurements in degrees from the magnetic compass
    Define_Constants;

    lambda_0 = GNSS_Solution.lambda_b(1);
    L_0 = GNSS_Solution.L_b(1);
    h = mean(GNSS_Solution.h_b); % check here for every GNSS data
    v_bar = (DR_Data(:,4)+DR_Data(:,5))/2 ;
    
    t = DR_Data(:,1);
    psi = psi_C; %DR_Data(:,7) * deg_to_rad;
    
    L_k(1) = L_0;
    lambda_k(1) = lambda_0;
    for k = 2:length(psi)
        v_N_bar(k) = 0.5*(cos(psi(k)) + cos(psi(k-1)))*v_bar(k);
        v_E_bar(k) = 0.5*(sin(psi(k)) + sin(psi(k-1)))*v_bar(k);
        
        [R_N,R_E] = Radii_of_curvature(L_k(k-1));
        L_k(k) = L_k(k-1) + (v_N_bar(k)*(t(k) - t(k-1))) / (R_N + h);
        lambda_k(k) = lambda_k(k-1) + (v_E_bar(k)*(t(k) - t(k-1)) / ((R_E + h) * cos(L_k(k)))); 
    end
    
    v_N(1) = v_bar(1) * cos(psi(1));
    v_E(1) = v_bar(1) * sin(psi(1));
    for k = 2:length(psi)
       v_N(k) = 1.7 * v_N_bar(k) - 0.7 * v_N(k-1);   
       v_E(k) = 1.7 * v_E_bar(k) - 0.7 * v_E(k-1);   
    end
    DR_Solution.L_b = L_k;
    DR_Solution.lambda_b = lambda_k;
    DR_Solution.Psi = psi;
    DR_Solution.v_E = v_E;
    DR_Solution.v_N = v_N;
    DR_Solution.t = DR_Data(:,1);

    figure
    plot(lambda_k*rad_to_deg, L_k*rad_to_deg, '-');
    title('Dead Reckoning Position Result');
    xlabel('Latitude (deg)');
    ylabel('Longitude (deg)');
    grid on;

    figure
    quiver(lambda_k*rad_to_deg ,L_k*rad_to_deg,v_E,v_N)
    grid on;
    xlabel('Longitude (deg)');
    ylabel('Latitude (deg)');
    title('Damped Instantaneous DR Velocity')
end