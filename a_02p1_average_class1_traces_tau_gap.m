clear
load 'all_sorted_traces_240614_to_240723_sorted_by_xb_6p5kbraid_v2_0p02turn.mat'
DAQ_all = a_DAQ_sorted_trace;

[row column] = size(DAQ_all);
num_of_traces = column/4;

selected_trace_num = [3:13];

selected_trace_count = 1;
for n = selected_trace_num
    length_of_data(selected_trace_count) = length(DAQ_all{:,4*n});
    
    turn_low_all_short(selected_trace_count) = DAQ_all{:,4*n-3}(1);
    turn_high_all_short(selected_trace_count) = DAQ_all{:,4*n-3}(end);
    
    
    figure(1)
    subplot(2,1,1)
    plot(DAQ_all{:,4*n-3},DAQ_all{:,4*n-2})%turn vs extension
    hold on
    subplot(2,1,2)
    plot(DAQ_all{:,4*n-1},DAQ_all{:,4*n})%turn vs torque
    hold on

    figure(9)
    plot(DAQ_all{:,4*n-3},DAQ_all{:,4*n-2})%turn vs extension
    hold on

    selected_trace_count = selected_trace_count + 1;
end

turn_low_cutoff = max(turn_low_all_short) + 1.8;
turn_high_cutoff = min(turn_high_all_short) - 1;
turn_with_smallest_range = DAQ_all{:,4*selected_trace_num(1)-3};
length_of_data_used = length(turn_with_smallest_range(turn_with_smallest_range > turn_low_cutoff & turn_with_smallest_range < turn_high_cutoff));
%use data start from 'turn_low_cutoff' and total length of 'length_of_data_used'
%%
selected_trace_count = 1;
for n = selected_trace_num
    turn_all_crop_short = DAQ_all{:,4*n-3}(DAQ_all{:,4*n-3} > turn_low_cutoff);
    extension_all_crop_short = DAQ_all{:,4*n-2}(DAQ_all{:,4*n-3} > turn_low_cutoff);
    
    turn_all_crop = DAQ_all{:,4*n-1}(DAQ_all{:,4*n-3} > turn_low_cutoff);
    torque_all_crop = DAQ_all{:,4*n}(DAQ_all{:,4*n-3} > turn_low_cutoff);
    
    turn_all_crop_short_final(:,selected_trace_count) = turn_all_crop_short(1:length_of_data_used);
    extension_all_crop_short_final(:,selected_trace_count) = extension_all_crop_short(1:length_of_data_used);
    turn_all_crop_final(:,selected_trace_count) = turn_all_crop(1:length_of_data_used);
    torque_all_crop_final(:,selected_trace_count) =  torque_all_crop(1:length_of_data_used);

    figure(2)
    subplot(2,1,1)
    plot(turn_all_crop_short_final(:,selected_trace_count),extension_all_crop_short_final(:,selected_trace_count))%turn vs extension
    hold on
    subplot(2,1,2)
    plot(turn_all_crop_final(:,selected_trace_count),torque_all_crop_final(:,selected_trace_count))%turn vs torque
    hold on
    
    %individual get tau_gap+ and tau_gap-
    turn_current = DAQ_all{:,4*n-1};
    torque_current = DAQ_all{:,4*n};
    figure(6)
    plot(turn_current(turn_current > 0.5 & turn_current < 10),torque_current(turn_current > 0.5 & turn_current < 10))
    hold on
    plot(turn_current(turn_current > -10 & turn_current < -0.5),torque_current(turn_current > -10 & turn_current < -0.5))
    hold on
    
    torque_gap_pos(selected_trace_count) = mean(torque_current(turn_current > 0.5 & turn_current < 10));
    torque_gap_neg(selected_trace_count) = mean(torque_current(turn_current > -10 & turn_current < -0.5));

    selected_trace_count = selected_trace_count + 1;
end

diff_torque_gap = torque_gap_pos - torque_gap_neg

[mean(diff_torque_gap) std(diff_torque_gap)]

