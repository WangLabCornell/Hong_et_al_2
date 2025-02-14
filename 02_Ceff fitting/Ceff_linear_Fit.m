function [fitresult, gof] = Ceff_linear_Fit(sigma_braid_for_fit, torque_braid_for_fit)
%CREATEFIT(SIGMA_BRAID_FOR_FIT,TORQUE_BRAID_FOR_FIT)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input: sigma_braid_for_fit
%      Y Output: torque_braid_for_fit
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 08-Aug-2024 20:37:37


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( sigma_braid_for_fit, torque_braid_for_fit );

% Set up fittype and options.
ft = fittype( 'k*x+b', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0.854099949273443 0.446026648055103];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
figure;
h = plot( fitresult, xData, yData );
legend( h, 'torque_braid_for_fit vs. sigma_braid_for_fit', 'untitled fit 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
xlabel( 'sigma_braid_for_fit', 'Interpreter', 'none' );
ylabel( 'torque_braid_for_fit', 'Interpreter', 'none' );
grid on
hold on


