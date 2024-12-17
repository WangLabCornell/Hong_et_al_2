1. plot_all_traces_setup4_braiding_v2_raw.m
Script to individually load the converted data including force, extension, turn, and torque for experiments described in Figure 2-4.

2. a_02_average_centered_hat_curves_3p0pN_for_braiding.m
subfunction: linearFit.m
Script to average extension vs. turn relation and torque vs. turn relation for all individual traces.

3. a_02_average_centered_hat_curves_3p0pN_for_braiding_v2.m
subfunction: Ceff_linear_Fit.m
Script to do a linear fit to individual traces to obtain the effective twist persistence length C_eff.

4. a_04_individual_hat_curves_3p0pN_for_braiding_torque.m
subfunction: linear_Fit.m
Script to obtain the post-buckling torque for individual traces.

5. a_05_buckling_transition_hat_curves_3p0pN_for_braiding.m
subfunction: buckling_transition_Fit_v2.m
Script to fit the buckling transition for individual traces.

6. a_06_individual_hat_curves_3p0pN_for_ext_slope.m
subfunction: ext_slope_Fit.m
Script to fit the post-buckling extension slope.

7. a_02p1_average_class1_traces_tau_gap.m
Script to obtain the torque gap for individual traces of the V substrate.

8. a_01_sort_traces_by_hat_tip_size_class1.m
subfunctions: Charvin_first_half_Fit.m ; torque_overshoot_Fit.m
Script to fit the hat tip size and the torque overshoot for individual traces of the U substrate.

9. supercoil_partition_2pN_gap_0pNnm_Ceff20nm_sta_at_0turn.m
Script to calculate supercoiling partition for the senerio of a/l = b/l = 0 (Figure 5a).

10. supercoiling_partition_2pN_tau_overshoot_significant.m
Script to calculate supercoiling partition for the senerio of a/l > 0 and b/l > 0 (Figure 5b).

11. supercoil_partition_2pN_gap_5pNnm.m
Script to calculate supercoiling partition for the senerio of a/l = 0 and b/l > 0 (Figure S5).