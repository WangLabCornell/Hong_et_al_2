load 'all_centered_traces_3pN_ext_no_smooth_tau_4turn.mat'
DAQ = a_centered_hat_curves;%turn extension torque



DAQ_all = [DAQ];

%%
[row column] = size(DAQ_all);
num_of_traces = column/4;

for n = 1:num_of_traces
    turn_all_for_ext = DAQ_all{:,4*n-3};
    torque_all_ext = DAQ_all{:,4*n-2};
    figure(2)
    plot(turn_all_for_ext,torque_all_ext)
    hold on
    
    turn_all_crop = DAQ_all{:,4*n-1};
    torque_all_crop = DAQ_all{:,4*n-0};
    
    
    turn_all_crop_ROI_center = turn_all_crop(turn_all_crop > -20 & turn_all_crop < 20);
    torque_all_crop_ROI_center = torque_all_crop(turn_all_crop > -20 & turn_all_crop < 20);
    
    [fitresult_tau_baseline, gof_tau_baseline] = linear_Fit(turn_all_crop_ROI_center, torque_all_crop_ROI_center)
    
    torque_all_crop_shifted = torque_all_crop - fitresult_tau_baseline.b;
    figure(1)
    plot(turn_all_crop,torque_all_crop_shifted)
    hold on


   
    a_turn_R_braid_buckle = turn_all_crop(turn_all_crop > 98);
    a_torque_R_braid_buckle = torque_all_crop_shifted(turn_all_crop > 98);

    figure(1)
    plot(a_turn_R_braid_buckle,a_torque_R_braid_buckle,'r')
    hold on

    buckling_torque(n) = mean(a_torque_R_braid_buckle)

end
%%

[mean(buckling_torque) std(buckling_torque)]