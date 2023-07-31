function outlier_labels = outlier_detection(detecting_epoch,H_G,delta_z)
    %%%%%%%%%%%%%%%%%%%%%%%
    % Input:  detecting_epoch, 1xN, current epoch of pseudo ranges
    %         H_G, Nx4, measurement matrix
    %         delta_z, Nx1, measurement innovation vector
    % Output: outlier_labels, nx2, index of outlier and residual vector
    %%%%%%%%%%%%%%%%%%%%%%%

    outlier_labels = [];
    sigma = 5; % measurement error standard deviation
    T = 6; % outlier detection threshold
    %% a) Compute the residuals vector
    residuals = (H_G/(H_G'*H_G)*H_G'-eye(size(detecting_epoch,2)-1))*delta_z;
    
    %% b) Compute the residuals convariance matrix
    C_v = (eye(size(detecting_epoch,2)-1)-H_G/(H_G'*H_G)*H_G')*sigma^2;

    %% c) Compute the normalised residuals and compare each with threshold
    for j = 1: length(residuals)
        if(abs(residuals(j)) > sqrt(C_v(j,j))*T)
            
            outlier_labels = [outlier_labels;[j,residuals(j)]];
        end
    end
end