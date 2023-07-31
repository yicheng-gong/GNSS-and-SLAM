function [x_est,P_matrix] = Initialise_GNSS_KF(pseudo_ranges_data,pseudo_range_rates_data )
    
    % compute least-square solution
    % and remove outlier
    while true
        [x_est,~,~,~,~,H_G,delta_z] = ...
            GNSS_LS_Solver(pseudo_ranges_data(1,:), ...
                           pseudo_ranges_data(2,:), ...
                           pseudo_range_rates_data(2,:), ...
                           zeros(8,1));
        outlier_labels = outlier_detection(pseudo_ranges_data(2,:),H_G,delta_z);
        if(~isempty(outlier_labels))
            [~,index] = max(abs(outlier_labels(:,2)));
            removel_index = outlier_labels(index,1);
            pseudo_ranges_data(:,removel_index+1) = [];
            pseudo_range_rates_data(:,removel_index+1) = [];
        else
            break
        end
    end
    
% Initialise error covariance matrix
P_matrix =  zeros(8);
P_matrix(1,1) = 100;
P_matrix(2,2) = 100;
P_matrix(3,3) = 100;
P_matrix(4,4) = 0.01;
P_matrix(5,5) = 0.01;
P_matrix(6,6) = 0.01;
P_matrix(7,7) = 100000^2;
P_matrix(8,8) = 200^2;

end