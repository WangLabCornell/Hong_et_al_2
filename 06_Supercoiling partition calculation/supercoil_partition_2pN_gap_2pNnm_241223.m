%single DNA
force_current = 2;%--------------------------------------------------
buckling_torque_single = 12.7825 * force_current.^0.632447;

load 'result_C_eff_vs_force_A43_C109_b6.0200.txt'
DAQ = result_C_eff_vs_force_A43_C109_b6_0200;
force_LUT = DAQ(:,1);
Ceff_LUT = DAQ(:,2);
Ceff_single_current = interp1(force_LUT,Ceff_LUT,force_current);
kBT = 4.09;
w0 = 2*pi/3.55;%3.55 nm = 0.338 nm/bp * 10.5 bp

sigma_s = buckling_torque_single/Ceff_single_current/kBT/w0;

sigma_upper = 0.12;%has to be smaller than sigma_p--------------------------------------------

%braided DNA
Ceff_braid = 25;%nm, treat as force-independent--------------------------------------------------
torque_gap = 2;%pNnm, has to be smaller than buckling_torque_single------------------------------

%regime 1
sigma1 = torque_gap / Ceff_single_current / kBT / w0;

%calculation

rho = 0.5;%replisome at midpoint
Lk0 = 14000/10.5;%full substrate, 14 kb
count = 1;
sigma_front_range = 0.00001:0.00001:sigma_upper;
for sigma_front = sigma_front_range
    if sigma_front <= sigma1%within torque gap region, assume no any turn
        torque_hinder_replisome(count) = Ceff_single_current * kBT * w0 * sigma_front;
        sigma_behind(count) = 0;
        turn_behind(count) = Ceff_single_current * kBT * w0 * sigma_front / (torque_gap/0.5);%< 0.5 turn, stiffness of braid is torque_gap/0.5turn
        turn_front(count) = sigma_front * (1 -rho) * Lk0;
    elseif sigma_front > sigma1 && sigma_front <= sigma_s
        torque_hinder_replisome(count) = Ceff_single_current * kBT * w0 * sigma_front;
        sigma_behind(count) = (Ceff_single_current * kBT * w0 * sigma_front - torque_gap) / Ceff_braid / kBT / w0;
        turn_behind(count) = sigma_behind(count) * rho * Lk0 + 0.5;%add the not-scaling 0.5 turn
        turn_front(count) = sigma_front * (1 -rho) * Lk0;
    elseif sigma_front > sigma_s
        torque_hinder_replisome(count) = buckling_torque_single;
        sigma_behind(count) = (buckling_torque_single - torque_gap) / Ceff_braid / kBT / w0;
        turn_behind(count) = sigma_behind(count) * rho * Lk0 + 0.5;%add the not-scaling 0.5 turn
        turn_front(count) = sigma_front * (1 -rho) * Lk0;
    end
    count = count + 1;
end


fraction_of_supercoiling_behind = turn_behind./(turn_behind + turn_front);
sigma_tot = (turn_behind + turn_front)/Lk0;

figure(1)
plot(sigma_tot,fraction_of_supercoiling_behind)
hold on
figure(2)
plot(sigma_tot,torque_hinder_replisome)
hold on
figure(3)
plot(sigma_front_range,torque_hinder_replisome,'r')
hold on

a_result = [sigma_tot' fraction_of_supercoiling_behind' torque_hinder_replisome'];