function Result_Plot(Method, L_b, lambda_b, color)
%% Positioning plot 
% figure
plot(lambda_b, L_b, '-', Color=color);
title(Method, ' Position');
xlabel('Latitude (deg)');
ylabel('Longitude (deg)');
grid on;

end
