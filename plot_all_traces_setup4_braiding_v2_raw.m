%Edit by 2023-02-10
%select the folder containing trace folders
clear
default_folder_path = '\\Yosemite\data\Setup4\Yifeng\240714_Yifeng_O template postbuckling v5 more turns_Processed'; % The defalt path for selecting data folders
select_folder_title = 'Please Select the Folder:'; % Title that displayed on the selecting folder window
current_folder=uigetdir(default_folder_path,select_folder_title); % current folder address
filePattern = dir(fullfile(current_folder, 'ProcessedTrace*.*')); % Find folder with folder name "ProcessedTrace"
[num_of_traces, sizey] = size(filePattern); % how many files (traces) are there
    
count_trace = 1;

x0_all = [0 0 0 0 0 0 0];
%centered by eye

for num = [2:6] %trace selection 2 3 4 5 6

    current_folder_sub = [char(filePattern(num).folder),'\',char(filePattern(num).name)];
    fileList = dir(fullfile(current_folder_sub,'00*.*'));% Returns a struct containing the matched name files

    
    %force clamp --------------------------------------------------------
        
        for num_of_force = 1
            if num_of_force == 1%0.5pN
                tasklist = 21:26;
            elseif num_of_force == 2%1pN
                tasklist = 32:37;
            elseif num_of_force == 3%2pN
                tasklist = 43:48;
            elseif num_of_force == 4%3pN
                tasklist = 54:59;
            end
            
            spin_count = 1;%a total of 6 spins
            sliding_window_short = 1;%0.0025s
            sliding_window = 801;%2s
            for kk = tasklist %
                cd(fileList(1).folder);%go to the current trace's path
                DAQ1 = importdata(fileList(kk).name);
                Data_all = DAQ1.data;
                
                cd(current_folder)
                Fz = Data_all(:,51);
                %Fz_smooth = movingmean(Fz,sliding_window,1,1);
                Zext = Data_all(:,52);
                %Zext_smooth = movingmean(Zext,sliding_window,1,1);
                Torque = Data_all(:,53);
                %Torque_smooth = movingmean(Torque,sliding_window,1,1);
                Turn = Data_all(:,54);
                %Turn_smooth  = movingmean(Turn,sliding_window,1,1);
                dzbead = Data_all(:,61);%negative values
                %dzbead_smooth = movingmean(dzbead,sliding_window,1,1);
                dz_cyl = -dzbead;
                
                h_trap = Zext - dzbead;
                %h_trap_smooth = Zext_smooth - dzbead_smooth;
                power = Data_all(:,41);
                %power_smooth = movingmean(power,sliding_window,1,1);
                
                
                %raw data
                figure(1)
                subplot(2,2,1)
                plot(Turn,Zext,'.')
                hold on
                figure(1)
                subplot(2,2,2)
                plot(Turn,Fz,'.')
                hold on
                figure(1)
                subplot(2,2,3)
                plot(Turn,Torque,'.')
                hold on
                figure(1)
                subplot(2,2,4)
                plot(Turn,power)
                hold on
                
%                 %smoothed data
%                 figure(2)
%                 subplot(2,1,1)
%                 plot(Turn_smooth,Zext_smooth)
%                 hold on
%                 subplot(2,1,2)
%                 plot(Turn_smooth,Torque_smooth)
%                 hold on
                
                %%
                %record the data of each wind
                Torque_all{spin_count} = Torque;
                Turn_all{spin_count} = Turn;
                Zext_all{spin_count} = Zext;
                Fz_all{spin_count} = Fz;
                h_trap_all{spin_count} = h_trap;
                dz_cyl_all{spin_count} = dz_cyl;
                
                spin_count = spin_count + 1;
                clear Turn
                clear Zext
                clear Fz
                clear Torque
                clear h_trap

            end
            %%
            %use wind 2 3 4 5 as continuous winding
            Turn_all{2} = flip(Turn_all{2});
            Turn_all{3} = flip(Turn_all{3});
            Turn_wind_2_3 = [Turn_all{3};Turn_all{2}];
            Turn_wind_4_5 = [Turn_all{4};Turn_all{5}];
            
            
            Zext_all{2} = flip(Zext_all{2});
            Zext_all{3} = flip(Zext_all{3});
            Zext_wind_2_3 = [Zext_all{3};Zext_all{2}];
            Zext_wind_4_5 = [Zext_all{4};Zext_all{5}];
            
            Torque_all{2} = flip(Torque_all{2});
            Torque_all{3} = flip(Torque_all{3});
            Torque_wind_2_3 = [Torque_all{3};Torque_all{2}];
            Torque_wind_4_5 = [Torque_all{4};Torque_all{5}];
            
            h_trap_all{2} = flip(h_trap_all{2});
            h_trap_all{3} = flip(h_trap_all{3});
            h_trap_wind_2_3 = [h_trap_all{3};h_trap_all{2}];
            h_trap_wind_4_5 = [h_trap_all{4};h_trap_all{5}];
            
            dz_cyl_all{2} = flip(dz_cyl_all{2});
            dz_cyl_all{3} = flip(dz_cyl_all{3});
            dz_cyl_wind_2_3 = [dz_cyl_all{3};dz_cyl_all{2}];
            dz_cyl_wind_4_5 = [dz_cyl_all{4};dz_cyl_all{5}];
            
            figure(2)
            subplot(2,2,1)
            plot(Turn_wind_2_3,Zext_wind_2_3)
            hold on
            plot(Turn_wind_4_5,Zext_wind_4_5)
            hold on
            figure(2)
            subplot(2,2,2)
            plot(Turn_wind_2_3,Torque_wind_2_3)
            hold on
            plot(Turn_wind_4_5,Torque_wind_4_5)
            hold on
            figure(2)
            subplot(2,2,3)
            plot(Turn_wind_2_3,h_trap_wind_2_3)
            hold on
            plot(Turn_wind_4_5,h_trap_wind_4_5)
            hold on
            figure(2)
            subplot(2,2,4)
            plot(Turn_wind_2_3,dz_cyl_wind_2_3)
            hold on
            plot(Turn_wind_4_5,dz_cyl_wind_4_5)
            hold on
            
            %%
            %smooth data
            Turn_wind_2_3_smooth_short = movingmean(Turn_wind_2_3,sliding_window_short,1,1);
            Turn_wind_4_5_smooth_short = movingmean(Turn_wind_4_5,sliding_window_short,1,1);
            Zext_wind_2_3_smooth_short = movingmean(Zext_wind_2_3,sliding_window_short,1,1);
            Zext_wind_4_5_smooth_short = movingmean(Zext_wind_4_5,sliding_window_short,1,1);
            
            Turn_wind_2_3_smooth = movingmean(Turn_wind_2_3,sliding_window,1,1);
            Turn_wind_4_5_smooth = movingmean(Turn_wind_4_5,sliding_window,1,1);
%             Zext_wind_2_3_smooth = movingmean(Zext_wind_2_3,sliding_window,1,1);
%             Zext_wind_4_5_smooth = movingmean(Zext_wind_4_5,sliding_window,1,1);
            Torque_wind_2_3_smooth = movingmean(Torque_wind_2_3,sliding_window,1,1);
            Torque_wind_4_5_smooth = movingmean(Torque_wind_4_5,sliding_window,1,1);
            
            h_trap_wind_2_3_smooth = movingmean(h_trap_wind_2_3,sliding_window,1,1);
            h_trap_wind_4_5_smooth = movingmean(h_trap_wind_4_5,sliding_window,1,1);
            
            dz_cyl_wind_2_3_smooth = movingmean(dz_cyl_wind_2_3,sliding_window,1,1);
            dz_cyl_wind_4_5_smooth = movingmean(dz_cyl_wind_4_5,sliding_window,1,1);
            
            figure(3)
            subplot(2,2,1)
            plot(Turn_wind_2_3_smooth_short,Zext_wind_2_3_smooth_short)
            hold on
            plot(Turn_wind_4_5_smooth_short,Zext_wind_4_5_smooth_short)
            hold on
            figure(3)
            subplot(2,2,2)
            plot(Turn_wind_2_3_smooth,Torque_wind_2_3_smooth)
            hold on
            plot(Turn_wind_4_5_smooth,Torque_wind_4_5_smooth)
            hold on
            figure(3)
            subplot(2,2,3)
            plot(Turn_wind_2_3_smooth,h_trap_wind_2_3_smooth)
            hold on
            plot(Turn_wind_4_5_smooth,h_trap_wind_4_5_smooth)
            hold on
            figure(3)
            subplot(2,2,4)
            plot(Turn_wind_2_3_smooth,dz_cyl_wind_2_3_smooth)
            hold on
            plot(Turn_wind_4_5_smooth,dz_cyl_wind_4_5_smooth)
            hold on

            figure(5)
            %subplot(2,2,1)
            plot(Turn_wind_2_3_smooth_short,Zext_wind_2_3_smooth_short)
            hold on
            plot(Turn_wind_4_5_smooth_short,Zext_wind_4_5_smooth_short)
            hold on
            xlim([-150 150])
            ylim([0 2500])
            
            %%
            %make sample size the same for each round - positive round 1 
            if length(Turn_wind_2_3_smooth) > length(Turn_wind_4_5_smooth)
                diff_error_1 = Turn_wind_2_3_smooth(1) - Turn_wind_4_5_smooth(1);
                diff_error_end = Turn_wind_2_3_smooth(end) - Turn_wind_4_5_smooth(end);
                delta_n = abs(length(Turn_wind_4_5_smooth) - length(Turn_wind_2_3_smooth));
                if abs(diff_error_1) <= abs(diff_error_end)
                    Turn_wind_2_3_smooth = Turn_wind_2_3_smooth(1:end-delta_n);
                    Torque_wind_2_3_smooth = Torque_wind_2_3_smooth(1:end-delta_n);
                    h_trap_wind_2_3_smooth = h_trap_wind_2_3_smooth(1:end-delta_n);
                    dz_cyl_wind_2_3_smooth = dz_cyl_wind_2_3_smooth(1:end-delta_n);
                    Turn_wind_2_3_smooth_short = Turn_wind_2_3_smooth_short(1:end-delta_n);
                    Zext_wind_2_3_smooth_short = Zext_wind_2_3_smooth_short(1:end-delta_n);
                elseif abs(diff_error_1) > abs(diff_error_end)
                    Turn_wind_2_3_smooth = Turn_wind_2_3_smooth(1+delta_n:end);
                    Torque_wind_2_3_smooth = Torque_wind_2_3_smooth(1+delta_n:end);
                    h_trap_wind_2_3_smooth = h_trap_wind_2_3_smooth(1+delta_n:end);
                    dz_cyl_wind_2_3_smooth = dz_cyl_wind_2_3_smooth(1+delta_n:end);
                    Turn_wind_2_3_smooth_short = Turn_wind_2_3_smooth_short(1+delta_n:end);
                    Zext_wind_2_3_smooth_short = Zext_wind_2_3_smooth_short(1+delta_n:end);
                end
            elseif length(Turn_wind_2_3_smooth) < length(Turn_wind_4_5_smooth)
                diff_error_1 = Turn_wind_2_3_smooth(1) - Turn_wind_4_5_smooth(1);
                diff_error_end = Turn_wind_2_3_smooth(end) - Turn_wind_4_5_smooth(end);
                delta_n = abs(length(Turn_wind_4_5_smooth) - length(Turn_wind_2_3_smooth));
                if abs(diff_error_1) <= abs(diff_error_end)
                    Turn_wind_4_5_smooth = Turn_wind_4_5_smooth(1:end-delta_n);
                    Torque_wind_4_5_smooth = Torque_wind_4_5_smooth(1:end-delta_n);
                    h_trap_wind_4_5_smooth = h_trap_wind_4_5_smooth(1:end-delta_n);
                    dz_cyl_wind_4_5_smooth = dz_cyl_wind_4_5_smooth(1:end-delta_n);
                    Turn_wind_4_5_smooth_short = Turn_wind_4_5_smooth_short(1:end-delta_n);
                    Zext_wind_4_5_smooth_short = Zext_wind_4_5_smooth_short(1:end-delta_n);
                elseif abs(diff_error_1) > abs(diff_error_end)
                    Turn_wind_4_5_smooth = Turn_wind_4_5_smooth(1+delta_n:end);
                    Torque_wind_4_5_smooth = Torque_wind_4_5_smooth(1+delta_n:end);
                    h_trap_wind_4_5_smooth = h_trap_wind_4_5_smooth(1+delta_n:end);
                    dz_cyl_wind_4_5_smooth = dz_cyl_wind_4_5_smooth(1+delta_n:end);
                    Turn_wind_4_5_smooth_short = Turn_wind_4_5_smooth_short(1+delta_n:end);
                    Zext_wind_4_5_smooth_short = Zext_wind_4_5_smooth_short(1+delta_n:end);
                end
            elseif length(Turn_wind_2_3_smooth) == length(Turn_wind_4_5_smooth)

            end
            
            %%
            %subtract viscous drag
            Turn_hat_curve = (Turn_wind_2_3_smooth + Turn_wind_4_5_smooth)/2;
            Torque_hat_curve = (Torque_wind_2_3_smooth + Torque_wind_4_5_smooth)/2;
            h_trap_hat_curve = (h_trap_wind_2_3_smooth + h_trap_wind_4_5_smooth)/2;
            dz_cyl_hat_curve = (dz_cyl_wind_2_3_smooth + dz_cyl_wind_4_5_smooth)/2;
            
            Turn_hat_curve_short = (Turn_wind_2_3_smooth_short + Turn_wind_4_5_smooth_short)/2;
            Zext_hat_curve_short = (Zext_wind_2_3_smooth_short + Zext_wind_4_5_smooth_short)/2;
            
%             figure(4)
%             subplot(2,1,1)
%             plot(Turn_hat_curve,Zext_hat_curve)
%             hold on
%             grid on
%             subplot(2,1,2)
%             plot(Turn_hat_curve,Torque_hat_curve)
%             hold on
%             grid on
            

            %%
            
            figure(4)
            subplot(2,1,1)
            plot(Turn_hat_curve_short-x0_all(num),Zext_hat_curve_short)
            hold on
            grid on

            figure(4)
            subplot(2,1,2)
            plot(Turn_hat_curve-x0_all(num),Torque_hat_curve)
            hold on
            grid on
            
            a_raw_DATA{:,count_trace*4-3} = Turn_hat_curve_short-x0_all(num);
            a_raw_DATA{:,count_trace*4-2} = Zext_hat_curve_short;
            a_raw_DATA{:,count_trace*4-1} = Turn_hat_curve-x0_all(num);
            a_raw_DATA{:,count_trace*4-0} = Torque_hat_curve;
            
            x0_all(num)
            



            count_trace = count_trace + 1;
    %force clamp end-------------------------------------------------------        
        end

end

