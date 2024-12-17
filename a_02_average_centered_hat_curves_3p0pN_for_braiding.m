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


for n = 1:num_of_traces
    figure(2)
    subplot(2,1,1)
    plot(DAQ_all_crop{:,4*n-3},DAQ_all_crop{:,4*n-2})
    hold on
    subplot(2,1,2)
    plot(DAQ_all_crop{:,4*n-1},DAQ_all_crop{:,4*n})
    hold on
    
    turn_all_crop_short(:,n) = DAQ_all_crop{:,4*n-3};
    extension_all_crop_short(:,n) = DAQ_all_crop{:,4*n-2};
    
    turn_all_crop(:,n) = DAQ_all_crop{:,4*n-1};
    torque_all_crop(:,n) = DAQ_all_crop{:,4*n};
end
%%

%average all traces
turn_average_short = mean(turn_all_crop_short');
extension_average_short = mean(extension_all_crop_short');

turn_average = mean(turn_all_crop');
torque_average = mean(torque_all_crop');

figure(3)
subplot(2,1,1)
plot(turn_average_short,extension_average_short)
hold on
subplot(2,1,2)
plot(turn_average,torque_average)
hold on


%%

%use 3 piece fitting, Figure(4)
turn_average_ROI = turn_average(turn_average > -2.5 & turn_average < 2.5);
torque_average_ROI = torque_average(turn_average > -2.5 & turn_average < 2.5);

[fitresult, gof] = linearFit(turn_average_ROI, torque_average_ROI);

torque_intecept = fitresult.b;

torque_average_shifted = torque_average - torque_intecept;

%%
figure(5)
subplot(2,1,1)
plot(turn_average_short,extension_average_short)
hold on
grid on
subplot(2,1,2)
plot(turn_average,torque_average_shifted)
hold on
grid on

a_final_hat_curve = {turn_average_short'  extension_average_short'  turn_average'  torque_average_shifted'};

