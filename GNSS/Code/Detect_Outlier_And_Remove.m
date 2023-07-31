function [satellite_name_list, pseudo_ranges, pseudo_range_rates ] = Detect_Outlier_And_Remove(x_pri,epoch, pseudo_ranges_data, pseudo_range_rates_data)
    while true
        % compute H matrix and z
        [~,~,~,~,~,H_G,delta_z] = ...
             GNSS_LS_Solver(pseudo_ranges_data(1,:), ...
                            pseudo_ranges_data(epoch,:), ...
                            pseudo_range_rates_data(epoch,:), ...
                            x_pri);
        
        % compute outlier index
        outlier_labels = outlier_detection(pseudo_ranges_data(epoch,:),H_G,delta_z);

        % check outlier: none, then stop
        if(isempty(outlier_labels))
            break
        end

        % check outlier: yes, remove the maximum
        [~,index] = max(abs(outlier_labels(:,2)));
        removel_index = outlier_labels(index,1);
        pseudo_ranges_data(:,removel_index+1) = [];
        pseudo_range_rates_data(:,removel_index+1) = [];
    end

    % store pseudo_ranges, pseudo_range_rates
    % and satellite name after outlier
    satellite_name_list = pseudo_ranges_data(1,2:end);
    pseudo_ranges = pseudo_ranges_data(epoch,2:end);
    pseudo_range_rates = pseudo_range_rates_data(epoch,2:end);
end