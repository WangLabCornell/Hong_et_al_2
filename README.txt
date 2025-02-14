- This repository contains a list of example scripts, based on MATLAB R2020a, for analysis on the traces measured by the angular optical trap (AOT). 

- MATLAB R2020a can be downloaded from MathWorks (https://www.mathworks.com/products/matlab.html) and installed with license. In our case, it was installed and run on a Window 10 operation system. MATLAB can be installed within 30 min. No special packages are required to run these scripts.

- To run these scripts, you may simply click run by opening the script with MATLAB, as the example dataset is placed within the same folder. Typically, the script can be run within seconds.

- Detailed function of each script:
01_Averaging traces measured by AOT
Data: 'all_centered_traces_3pN_ext_no_smooth_tau_4turn.mat'
Description of data: Individual traces of torsional measurements for the 'O substrate' held at 3 pN.
Main analysis script: 'a_02_average_centered_hat_curves_3p0pN_for_braiding.m'
Description of the main script: Averaging all the individual traces to obtain the extension vs. turn relation and torque vs. turn relation while removing the torque baseline. The output refers to one of the trace shown in Fig. 2b. Similar analysis applys to Fig. 3a, Fig. 3b, and Fig. 4c.

02_Ceff fitting
Data: 'all_centered_traces_3pN_ext_no_smooth_tau_4turn.mat'
Description of data: Individual traces of torsional measurements for the 'O substrate' held at 3 pN.
Main analysis script: 'a_02_average_centered_hat_curves_3p0pN_for_braiding_v2.m'
Description of the main script: Performing a linear fitting for the effective twist persistence length C_eff to individual traces. The output (mean value ± s.d.) refers to one of the data point shown in Fig. 2c. Similar analysis applys to other data points in Fig. 3c.

03_Obtaining buckling torque of a DNA braid
Data: 'all_centered_traces_3pN_ext_no_smooth_tau_4turn.mat'
Description of data: Individual traces of torsional measurements for the 'O substrate' held at 3 pN.
Main analysis script: 'a_04_individual_hat_curves_3p0pN_for_braiding_torque.m'
Description of the main script: Obtaining the post-buckling torque of a DNA braid via averaging the torque value measured after the DNA braid buckles. The output (mean value ± s.d.) refers to one of the data point shown in Fig. 2d (top panel). Similar analysis applys to other data points in Fig. 3d (top panel).

04_Fitting buckling transition of a DNA braid
Data: 'all_centered_traces_3pN_ext_no_smooth_tau_4turn.mat'
Description of data: Individual traces of torsional measurements for the 'O substrate' held at 3 pN.
Main analysis script: 'a_05_buckling_transition_hat_curves_3p0pN_for_braiding.m'
Description of the main script: Obtaining the buckling transition position by fitting the individual extension-turn relations with a quadratic funtion. The output (mean value ± s.d.) refers to one of the data point shown in Fig. 2d (bottom panel). Similar analysis applys to other data points in Fig. 3d (bottom panel).

05_Fitting post buckling extension slope
Data: 'all_centered_traces_3pN_ext_no_smooth_tau_4turn.mat'
Description of data: Individual traces of torsional measurements for the 'O substrate' held at 3 pN.
Main analysis script: 'a_06_individual_hat_curves_3p0pN_for_ext_slope.m'
Description of the main script: Obtaining the post-buckling extension by linearly fitting the individual extension-turn relations after the DNA braid buckles. The output (mean value ± s.d.) refers to one of the data point shown in Fig. 2e.

06_Supercoiling partition calculation
Data: 'result_C_eff_vs_force_A43_C109_b6.0200.txt'
Description of data: Simulation of the force vs. C_eff relation for a single DNA molecule based on the Bouchiat and M´ezard (BM) theory (X. Gao, et al., Torsional Stiffness of Extended and Plectonemic DNA, Physical Review Letters, 2021).
Main analysis script: 'supercoil_partition_2pN_gap_0pNnm_Ceff20nm_sta_at_0turn.m' for supercoiling partitioning calculation for a/l = b/l = 0; 'supercoiling_partition_2pN_tau_overshoot_significant_241223.m'  for supercoiling partitioning calculation for a/l >> 0 ,b/l >> 0; 'supercoil_partition_2pN_gap_2pNnm_241223.m' for supercoiling partitioning calculation for a/l = 0, b/l >> 0.
Description of the main script: Calcuating the geometry-dependent supercoiling partitioning assuming the replisome translocates to the midpoint of a 14 kb DNA substrate. Each output refers to Fig. 5a, Fig. 5b, and Fig. S5.
