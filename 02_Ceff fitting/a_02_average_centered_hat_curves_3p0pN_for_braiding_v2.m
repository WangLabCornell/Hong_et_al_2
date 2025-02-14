load 'all_centered_traces_3pN_ext_no_smooth_tau_4turn.mat'
DAQ = a_centered_hat_curves;%turn extension torque



DAQ_all = [DAQ];

%%
%find the length used for averaging
[row column] = size(DAQ_all);
num_of_traces = column/4;

for n = 1:num_of_traces
    length_of_data(n) = length(DAQ_all{:,4*n});
    
    turn_low_all_short(n) = DAQ_all{:,4*n-3}(1);
    turn_high_all_short(n) = DAQ_all{:,4*n-3}(end);
    
    turn_low_all(n) = DAQ_all{:,4*n-1}(1);
    turn_high_all(n) = DAQ_all{:,4*n-1}(end);
    
    figure(1)
    subplot(2,1,1)
    plot(DAQ_all{:,4*n-3},DAQ_all{:,4*n-2})%turn vs extension
    hold on
    subplot(2,1,2)
    plot(DAQ_all{:,4*n-1},DAQ_all{:,4*n})%turn vs torque
    hold on
end

[min_length index0] = min(length_of_data);
length_used = length_of_data(index0);%use the minimum length as the length of each trace

[max_turn_low_short index1] = max(turn_low_all_short);
[max_turn_low index2] = max(turn_low_all);

%make each trace have the same length
for n = 1:num_of_traces
    [minValue_short closest_index_short] = min(abs(DAQ_all{:,4*n-3} - max_turn_low_short));
    [minValue closest_index] = min(abs(DAQ_all{:,4*n-1} - max_turn_low));
    
    closest_index_all_short(n) = closest_index_short;%for extension
    closest_index_all(n) = closest_index;%for torque
end

%%
%torque data

length_for_average_short = min_length - (max(closest_index_all_short)-1);%this is the length to crop to for extension
length_for_average = min_length - (max(closest_index_all)-1);%this is the length to crop to for torque

for n = 1:num_of_traces
    DAQ_all_crop{:,4*n-3} = DAQ_all{:,4*n-3}(closest_index_all_short(n):(closest_index_all_short(n)+length_for_average_short-1));
    DAQ_all_crop{:,4*n-2} = DAQ_all{:,4*n-2}(closest_index_all_short(n):(closest_index_all_short(n)+length_for_average_short-1));
    
    DAQ_all_crop{:,4*n-1} = DAQ_all{:,4*n-1}(closest_index_all(n):(closest_index_all(n)+length_for_average-1));
    DAQ_all_crop{:,4*n} = DAQ_all{:,4*n}(closest_index_all(n):(closest_index_all(n)+length_for_average-1));
end
%%

count_trace = 1;
for n = 1:num_of_traces
    
    turn_average = DAQ_all_crop{:,4*n-1};
    torque_average_unshifted = DAQ_all_crop{:,4*n};

    
    %%
    %fit Rbraid
    Lk0 = 7092/10.5;%------------------------------------------------------
    kBT = 4.09;%pNnm
    w0 = 2*pi/3.55;%3.55 nm = 0.338nm/bp*10.5bp
    sigma_braid = turn_average/Lk0;
    figure(6)
    plot(sigma_braid,torque_average_unshifted)
    hold on
    
    
    
    upper_limit_for_fit = 0.05;%--------------------------------------------
    sigma_braid_for_fit = sigma_braid(sigma_braid > 5/Lk0 & sigma_braid < upper_limit_for_fit);
    torque_braid_for_fit = torque_average_unshifted(sigma_braid > 5/Lk0 & sigma_braid < upper_limit_for_fit);
    a_braiding_torque_fitting_range = [sigma_braid_for_fit' torque_braid_for_fit'];
    
    [fitresult_Ceff_braid, gof_Ceff_braid] = Ceff_linear_Fit(sigma_braid_for_fit, torque_braid_for_fit);
    torque_slope_in_sigma = fitresult_Ceff_braid.k;
    
    Ceff_braid_fit(count_trace) = torque_slope_in_sigma/kBT/w0
    count_trace = count_trace + 1;
end
%%
Ceff_Rbraid_mean = mean(Ceff_braid_fit);
Ceff_Rbraid_std = std(Ceff_braid_fit);
a_Ceff_final = [Ceff_Rbraid_mean Ceff_Rbraid_std]