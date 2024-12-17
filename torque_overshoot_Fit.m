function [fitresult, gof] = torque_overshoot_Fit(turn_average_ROI, torque_average_ROI)
%CREATEFIT(TURN_AVERAGE_ROI,TORQUE_AVERAGE_ROI)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input: turn_average_ROI
%      Y Output: torque_average_ROI
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 12-Jul-2024 14:25:38


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( turn_average_ROI, torque_average_ROI );

% Set up fittype and options.
ft = fittype( 'A*sin(2*pi*x)+C', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [0 -Inf];
opts.StartPoint = [0.823457828327293 0.694828622975817];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
figure;
h = plot( fitresult, xData, yData );
legend( h, 'torque_average_ROI vs. turn_average_ROI', 'untitled fit 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
xlabel( 'turn_average_ROI', 'Interpreter', 'none' );
ylabel( 'torque_average_ROI', 'Interpreter', 'none' );
grid on
hold on

