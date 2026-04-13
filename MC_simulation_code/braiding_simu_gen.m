function [vj1, vj2]= braiding_simu_gen(n0, F, Lp, dtopbp,dbotbp, dLk0, d0, nc, trials, kmax, ftrial, v1_initial, v2_initial, seg, conc, kt, ko,kd,repeat)
% edited from braiding_simu_otemp_cont_conc_gen2_kt

rcutoff  = ftrial.GridVectors{1}(end);

toln = 1e-10;
kT=4.1; % pN nm 
rise = 0.338;       % nm
L0 = n0*rise;       % contour length in nm for both strands
n = round(L0/seg);  % number of segments, for both strand, n is even
if mod(n,2) ~= 0
    n = n+1;      % if n is odd increase length by 1 segment of 10 nm
end

Zfar = 200;        % distance to infinity, in segment

% %% distance at bottom, 2 chains are anchored along x, distance by 2*e
% %% bottom anchoring of chain 2 is on +x compared to bottom anchoring of chain 1
et = (dtopbp*0.338)/2; % 'e' in nm
et = et/seg; % e in seg units
dt = 2*et; % top distance

eb = (dbotbp*0.338)/2; % 'e' in nm
eb = eb/seg; % e in seg units
db = 2*eb; % bottom distance

% spring constant for additional energy terms
fact=Lp/(2*seg);
k_b = kT * fact;     % homogenous bending stiffness
% ko = kT * 20;     % to align the chain to +z
% kd = kT * 20;     % to fix the distance between 2 top anchoring points
% kt = kT * 20;     % to align the direction of the two connecting vectors at bottom and at top

vj1= v1_initial;
vj2 = v2_initial;
dv_extend1 = [ -Zfar, vj1(2,end), vj1(3,end) ;-Zfar,0,0]';   % wings of a square
dv_extend2 = [  Zfar, vj2(2,end), vj2(3,end) ; Zfar,0,0]';   % wings of a square

v1_extend = [vj1 , dv_extend1 , vj1(:,1)];
v2_extend = [vj2 , dv_extend2 , vj2(:,1)];
v1_extendxz = [v1_extend(3,:); v1_extend(1,:); v1_extend(2,:)] ;
v2_extendxz = [v2_extend(3,:); v2_extend(1,:); v2_extend(2,:)] ;    
choices = nchoosek(1:length(v1_extendxz)+length(v2_extendxz)-1,2);
%choosing_r2 = [choices;fliplr(choices)];
lengthall = length(vj1)+length(vj2)-1;
choices_rep = nchoosek(1:length(vj1)+length(vj2)-1,2);
iarr = 2:lengthall;
jarr = 1:lengthall-1;
test1 = (iarr-1).*((lengthall+1)-iarr/2)+jarr-iarr+2;
test1 = [1,test1(1:end-1)];
test0 = 0:lengthall-2;
test1 = test1 - test0;
%choosearr = nchoosek(1:length(r)-1,2);
choices_rep(test1,:) =[];
index2 = choices_rep(:,1)==length(vj1) | choices_rep(:,2)==length(vj1);
choices_rep(index2,:) =[];


dLk = Alex_link_speedx1v3(v1_extendxz, v2_extendxz, -1.1 + toln, -1.1 + toln,choices);
disp(['GIVEN : ' num2str(dLk0)])
disp(['Del (t=-1.1) : ' num2str(dLk)])
disp(['Extension : ' num2str(mean([vj1(3,end) vj2(3,end)]))])
disp(['Hard Wall Energy vector (d0=' num2str(d0) ') : ' num2str(E_coulomb4_vector(vj1,vj2,5,seg))])
disp(['Repulsion Energy vector (d0=' num2str(d0) ') : ' num2str(E_debye_huckel_no_neighbor_cubic(vj1,vj2,seg, rcutoff, ftrial,choices_rep,index2))])


hard_core_count1 = sum(vj1(3,:)>=0 & vj1(3,:) <= vj1(3,end));
hard_core_count2 = sum(vj2(3,:)>=0 & vj2(3,:) <= vj2(3,end)); 
if hard_core_count1 ~= n+1 || hard_core_count2 ~= n+1 % energy accepted    
    disp('flag')
end

vj1 = v1_initial;
vj2 = v2_initial;

% open binary files to save data and configurations
% fidbin1 = fopen(['intermediate_data_bin_' num2str(n0) 'bp_seg_' num2str(seg) 'nm_lp_' num2str(Lp) 'nm_kt_' num2str(kt) 'pNnm_ko_' num2str(ko) 'pNnm_kd_' num2str(kd) 'pNnm_dt_' num2str(dt) 'nmdb_' num2str(db) 'nm_' num2str(F) '_pN_Na' num2str(conc) 'mM_repeat' num2str(repeat) 'Cat' num2str(dLk0) 'turns.bin' ],'a');         % 'a' = append data
% fidbinc1 = fopen(['config1_data_bin_' num2str(n0) 'bp_seg_' num2str(seg) 'nm_lp_' num2str(Lp) 'nm_kt_' num2str(kt) 'pNnm_ko_' num2str(ko) 'pNnm_kd_' num2str(kd) 'pNnm_dt_' num2str(dt) 'nmdb_' num2str(db) 'nm_' num2str(F) '_pN_Na' num2str(conc) 'mM_repeat' num2str(repeat) 'Cat' num2str(dLk0) 'turns.bin' ],'a');         % 'a' = append data
% fidbinc2 = fopen(['config2_data_bin_' num2str(n0) 'bp_seg_' num2str(seg) 'nm_lp_' num2str(Lp) 'nm_kt_' num2str(kt) 'pNnm_ko_' num2str(ko) 'pNnm_kd_' num2str(kd) 'pNnm_dt_' num2str(dt) 'nmdb_' num2str(db) 'nm_' num2str(F) '_pN_Na' num2str(conc) 'mM_repeat' num2str(repeat) 'Cat' num2str(dLk0) 'turns.bin' ],'a');         % 'a' = append data

    %% Ei
    tau2 = vj2(:,end) - vj1(:,end);
    tau2_n = sqrt(tau2(1)^2 + tau2(2)^2 + tau2(3)^2);
    tau2_n_2d = sqrt(tau2(1)^2 + tau2(2)^2);
    beta = acos(tau2(1)/tau2_n_2d);
    mod_beta = sign(tau2(2))*beta;
    Ei=Ebend_vector(vj1,k_b) + Ebend_vector(vj2,k_b) - F*0.5*(vj1(3,end)+vj2(3,end))*seg + ...
                0.5*ko * tau2(3)^2 + 0.5*kd* (tau2_n-dt)^2 * seg^2 + 0.5 * kt * mod_beta^2 + E_debye_huckel_no_neighbor_cubic(vj1,vj2,seg, rcutoff, ftrial,choices_rep,index2);

    k = 0;
    fail_cat =0;
    fail_rep =0;
    fail_ener =0;
    
    for trial =1 :trials
    %% for pivot
    theta_stretch= -0.5236 + 1.0472*rand(nc,1); % +/-30 deg
    pp_stretch=randi([1,n],nc,1);
    di_stretch=randn(nc,3);         %random direction
    di_stretch=di_stretch./repmat(sqrt(sum(di_stretch.*di_stretch,2)),[1 3]);
    rm_stretch=rand(1,nc);

    %% for CR
    phi= -0.8727 + 1.7453*rand(nc,1); %random rotation angles between -50 and +50 degrees
    pp2=randi([1,n+1],round(1.2*nc),2);
    pps1 = sort(pp2,2)';
    pps1 = pps1(:,diff(pps1)>=3);
    pps1 = pps1(:,1:nc)';
    r_choice=rand(2,nc); 
    
    %% for reptation, only do reptation for over 3 segments
    pp2=randi([2,n],round(20*nc),2);
    ppr1 = sort(pp2,2)';
    ppr1 = ppr1(:,diff(ppr1)>=3 & diff(ppr1)<=10);
    ppr1 = ppr1(:,1:nc)';
    r_direction=rand(1,nc);
    
    index_arr = false(nc,1);
    pivot_arr =false(nc,1);
    cr_arr =false(nc,1);
    rep_arr =false(nc,1);

    for j=1:nc
        v1=vj1; v2 = vj2;
        
        if r_choice(1,j)< 0.5
            if r_choice(2,j)< 0.5
                v1(:,pp_stretch(j)+1:end)=pivot_single(v1,pp_stretch(j),di_stretch(j,:),theta_stretch(j)); % pivot
                pivot_arr(j)= true;
            elseif r_choice(2,j)>= 0.5 && r_choice(2,j)< 0.75
                v1(:,pps1(j,1)+1:pps1(j,2)-1) = CR(v1,pps1(j,1),pps1(j,2),phi(j)); % crankshaft
                cr_arr(j) =true;
            else
                v1 = reptation(vj1,ppr1(j,1),ppr1(j,2),r_direction(j));
                rep_arr(j) =true;
            end
        else
            if r_choice(2,j)< 0.5
                v2(:,pp_stretch(j)+1:end)=pivot_single(v2,pp_stretch(j),di_stretch(j,:),theta_stretch(j)); % pivot
                pivot_arr(j)= true;
            elseif r_choice(2,j)>= 0.5 && r_choice(2,j)< 0.75
                v2(:,pps1(j,1)+1:pps1(j,2)-1) = CR(v2,pps1(j,1),pps1(j,2),phi(j)); % crankshaft
                cr_arr(j) =true;            
            else 
                v2 = reptation(vj2,ppr1(j,1),ppr1(j,2),r_direction(j));
                rep_arr(j) =true;
            end 
        end
        
            tau2 = v2(:,end) - v1(:,end);
            %tau2(3) = 0;
            tau2_n = sqrt(tau2(1)^2 + tau2(2)^2 + tau2(3)^2);
            tau2_n_2d = sqrt(tau2(1)^2 + tau2(2)^2);
            beta = acos(tau2(1)/tau2_n_2d);
            mod_beta = sign(tau2(2))*beta;
            
            Ef=Ebend_vector(v1,k_b) + Ebend_vector(v2,k_b) - F*0.5*(v1(3,end)+v2(3,end))*seg + ...
                        0.5*ko * tau2(3)^2 + 0.5*kd* (tau2_n-dt)^2 * seg^2 + 0.5 * kt * mod_beta^2 + E_debye_huckel_no_neighbor_cubic(v1,v2,seg, rcutoff, ftrial, choices_rep,index2);
            
            if abs(Ef) > 10e5
               fail_rep =  fail_rep+1;
            elseif rm_stretch(j) > exp(-(Ef-Ei)/kT) 
               fail_ener =  fail_ener+1;
            end
            
            hard_core_count1 = sum(v1(3,:)>=0 & v1(3,:) <= v1(3,end));
            hard_core_count2 = sum(v2(3,:)>=0 & v2(3,:) <= v2(3,end)); 

            %Metropolis algorithm
            if rm_stretch(j) < exp(-(Ef-Ei)/kT) && hard_core_count1 == n+1 && hard_core_count2 == n+1 % energy accepted    
                    %% checking linking number
                    % closing the loops
                    dv_extend1 = [ -Zfar, v1(2,end), v1(3,end) ;-Zfar,0,0]'; % wings of a square
                    dv_extend2 = [Zfar, v2(2,end), v2(3,end) ;Zfar,0,0]';   % wing of a square
                    v1_extend = [v1 , dv_extend1 , v1(:,1)];
                    v2_extend = [v2 , dv_extend2 , v2(:,1)];
                    v1_extendxzn = [v1_extend(1,:); v1_extend(3,:); v1_extend(2,:)] ;
                    v2_extendxzn = [v2_extend(1,:); v2_extend(3,:); v2_extend(2,:)] ;    
                    % calculating the linking number    
                    try
                         Cf = Alex_link_speedx1v3(v1_extendxzn, v2_extendxzn, -1.1 + toln,-1.1 + toln, choices);
                    catch
                        warning('Problem using function.  Assigning a value of dLk-1.');
                        Cf = dLk0-1;
                    end
                    
                    % selecting only conformation with correct linking number  
                    if  mod(log(Cf/dLk)/log(1.1),1)<0.00005 ||   (1- mod(log(Cf/dLk)/log(1.1),1))<0.00005 
                        index_arr(j) =true;
                        k=k+1;  
                        vj1=v1;
                        vj2=v2;
                        Ei=Ef;

                        if rem(k,2000) == 0
                          disp('Progress');
                          %disp(k);
                          disp(j);
                          dLk1 = Alex_link_speedx1v3(v1_extendxzn, v2_extendxzn, -1 + toln, -1 + toln , choices);
                           disp((vj1(3,end)+vj2(3,end))/2);

                            plot3(vj1(1,:),vj1(2,:),vj1(3,:), 'r')
                            hold on
                            plot3(vj2(1,:),vj2(2,:),vj2(3,:), 'b')
                            hold off
                            axis equal 
                            drawnow;

                          if abs(dLk1 -dLk0)<1e-2
                            % fwrite(fidbinc1,vj1,'double');
                            % fwrite(fidbinc2,vj2,'double');
                            interdata = [vj1(3,end) vj2(3,end) mod_beta dLk1 Cf];
                            % fwrite(fidbin1,interdata,'double');
                          end
                        end
                    else
                        fail_cat =fail_cat+1;
                    end
            end
            if k>kmax
                break
            end
    end
            if k>kmax
                break
            end
  %  toc
    end

end

function detAlex = Alex_link_speedx1v3(r1, r2, s, t, choosing_r)
% calculation of alexander polynomial for 2 linked trajectories r1 and r2
% r1 and r2 are two closed loops, N1, and N2 are the number of segments, 
% Algorithm for Alexander_polynomial_calculation_linking_2 from Vologodskii and Rybenkov, Phys Chem Chem Phys, 2009
% algorithm 1 and 2 have similar speed
% The linking number is evaluated at s = 1 + tol and t = 1 + tol


tol = 1e-10;

N1 = size(r1,2)-1; % this is the first closed loop, N1 is the number of segment, N1+1 is the total number of points, r1_(N1+1) == r_1(1)
N2 = size(r2,2)-1; % this is the second closed loop, N2 is the number of segment,N2+1 is the total number of points, r2_(N2+1) == r_2(1)

r_combined = [r1(:,1:N1+1), r2(:,1:N2+1)];
d = diff(r_combined')';

aa = r_combined(:,choosing_r(:,1));
bb = r_combined(:,choosing_r(:,2));  
cc = r_combined(:,choosing_r(:,1)+1);  
dd = r_combined(:,choosing_r(:,2)+1);  


o1a = sign(((cc(2,:) - aa(2,:)) .* (bb(1,:) - cc(1,:))) - ((cc(1,:) - aa(1,:)) .* (bb(2,:) - cc(2,:)))); 
o2a = sign(((cc(2,:) - aa(2,:)) .* (dd(1,:) - cc(1,:))) - ((cc(1,:) - aa(1,:)) .* (dd(2,:) - cc(2,:)))); 
o3a = sign(((dd(2,:) - bb(2,:)) .* (aa(1,:) - dd(1,:))) - ((dd(1,:) - bb(1,:)) .* (aa(2,:) - dd(2,:)))); 
o4a = sign(((dd(2,:) - bb(2,:)) .* (cc(1,:) - dd(1,:))) - ((dd(1,:) - bb(1,:)) .* (cc(2,:) - dd(2,:)))); 

cond0 = choosing_r(:,1)== N1+1;
cond1 = abs(choosing_r(:,1) - choosing_r(:,2))<=1;
cond2 = choosing_r(:,1)== N1 & choosing_r(:,2)== 1;
cond3 = choosing_r(:,1)== 1 & choosing_r(:,2)== N1;
cond4 = choosing_r(:,1)== N1+N2+1 & choosing_r(:,2)== N1+2;
cond5 = choosing_r(:,1)== N1+2 & choosing_r(:,2)== N1+N2+1;
cond6 = choosing_r(:,2)== N1+1 ;

comb_conda = cond0 | cond1 | cond2 | cond3 | cond4 | cond5 | cond6; 
comb_conda = (~comb_conda);

intersecting = o1a(comb_conda) ~= o2a(comb_conda) & o3a(comb_conda) ~= o4a(comb_conda);
valid_choice1 = choosing_r(comb_conda,:);
valid_choice2 = valid_choice1(intersecting,:);

valid_choice2 =[valid_choice2;fliplr(valid_choice2)];

r_arr1_int(:,:,1) = r_combined(:,valid_choice2(:,1)) ;
r_arr1_int(:,:,2) = r_combined(:,valid_choice2(:,2)) ;

r_arr2_int(:,:,1) = r_combined(:,valid_choice2(:,1)+1) ;
r_arr2_int(:,:,2) = r_combined(:,valid_choice2(:,2)+1) ;

A1 = r_arr2_int(2,:,1)-r_arr1_int(2,:,1);
B1 = r_arr1_int(1,:,1)-r_arr2_int(1,:,1);
C1 = A1.*r_arr1_int(1,:,1)+ B1.*r_arr1_int(2,:,1);

A2 = r_arr2_int(2,:,2)-r_arr1_int(2,:,2);
B2 = r_arr1_int(1,:,2)-r_arr2_int(1,:,2);
C2 = A2.*r_arr1_int(1,:,2)+B2.*r_arr1_int(2,:,2);

detr = A1.*B2 - A2.*B1;
pt(1,:) = (B2.*C1 - B1.*C2)./detr;
pt(2,:) = (A1.*C2 - A2.*C1)./detr;

fd_r2 = sqrt( (pt(1,:)-r_arr1_int(1,:,2)).^2+(pt(2,:)-r_arr1_int(2,:,2)).^2 )./sqrt( (r_arr2_int(1,:,2)-r_arr1_int(1,:,2)).^2+(r_arr2_int(2,:,2)-r_arr1_int(2,:,2)).^2 );
fd_r1 = sqrt( (pt(1,:)-r_arr1_int(1,:,1)).^2+(pt(2,:)-r_arr1_int(2,:,1)).^2 )./sqrt( (r_arr2_int(1,:,1)-r_arr1_int(1,:,1)).^2+(r_arr2_int(2,:,1)-r_arr1_int(2,:,1)).^2 );

zl1 = ((pt(1,:)-r_arr1_int(1,:,1))./(r_arr2_int(1,:,1)-r_arr1_int(1,:,1))).*(r_arr2_int(3,:,1)-r_arr1_int(3,:,1)) + r_arr1_int(3,:,1);
zl2 = ((pt(1,:)-r_arr1_int(1,:,2))./(r_arr2_int(1,:,2)-r_arr1_int(1,:,2))).*(r_arr2_int(3,:,2)-r_arr1_int(3,:,2)) + r_arr1_int(3,:,2);

underpass_ind = zl1<zl2;

j_countu = sum(underpass_ind);
distance_of_under_u = fd_r1(underpass_ind) + valid_choice2(underpass_ind,1)';
distance_of_over_u = fd_r2(underpass_ind) + valid_choice2(underpass_ind,2)';
z_cross = d(1,valid_choice2(underpass_ind,2)).*d(2,valid_choice2(underpass_ind,1)) - d(2,valid_choice2(underpass_ind,2)).*d(1,valid_choice2(underpass_ind,1));  
type_under = sign(z_cross);


% if distance_of_under_u or distance_of_over_u is smaller than N1+1 then the crossing belongs to r1, otherwise it belongs to r
if j_countu >= 1
    % remove unnecessary elements
    distance_of_over_u(j_countu+1:end) = []; 
    distance_of_under_u(j_countu+1:end) = []; 
    type_under(j_countu+1:end) = [];
    
    M = sum(distance_of_under_u <= N1+1); % number of crossing in the contour r1
    N = j_countu; % total number of crossing

    if M == 0 || N == M
        detAlex = 0;
    else
        index_under_all = 1 : j_countu;
        [~,index_under] = sort(distance_of_under_u);
        distance_of_over_u = distance_of_over_u(index_under);
        distance_of_under_u = distance_of_under_u(index_under);
        type_under = type_under(index_under);

        % assinging generator
        generator = zeros(N,2); % the first collumn is the min, the second collumn is the max except for the first geenrator
        for i = 1 : M
            i2 = i+1;
            if i2 > M
                i2 = 1;
            end
            if abs(distance_of_under_u(i) - distance_of_under_u(i2)) < tol
                generator(i2,:)= [ 1 N1+1];
            else
                generator(i2,:)= [ distance_of_under_u(i) distance_of_under_u(i2)];
            end
        end
        
        if N >= M+1
            for i = M+1 : N
                i2 = i+1;
                if i2 > N
                    i2 = M+1;
                end
                if abs(distance_of_under_u(i) - distance_of_under_u(i2)) < tol
                    generator(i2,:)= [ N1+2 N1+N2+2];
                else
                    generator(i2,:)= [ distance_of_under_u(i) distance_of_under_u(i2)];
                end
            end
        end

        % finding index of the overpass according to the underpass
        index_generator = zeros(1,N);
        for i = 1 : N
            for j = 1 : N
                if generator(j,2) < generator(j,1)
                    if generator(j,1) >= N1+2
                        if (distance_of_over_u(i) > generator(j,1) && distance_of_over_u(i) <= N1+1 + N2 + 1 + tol)  || (distance_of_over_u(i) < generator(j,2) && distance_of_over_u(i) >= N1 + 2) 
                            index_generator(i) = j;
                        end 
                    else
                        if (distance_of_over_u(i) > generator(j,1) && distance_of_over_u(i) <= N1+1 + tol)  || (distance_of_over_u(i) < generator(j,2) && distance_of_over_u(i) >= 1)
                            index_generator(i) = j;
                        end
                    end
                else
                    if distance_of_over_u(i) > generator(j,1)  && distance_of_over_u(i) < generator(j,2)
                        index_generator(i) = j;
                    end
                end
            end
        end

        Ma_t = zeros(N,N);  
        %s = -1; t = -1;
        %s = 1+tol; t = 1+tol;
        for i = 1 : N
            k = index_under_all(i);
            ii = index_generator(i);
            if  k == N
                    kp1 = M+1;
            elseif k == M
                    kp1 = 1;
            else
                    kp1 = k+1;
            end

            if k <= M && M >1
                if ii == k || ii == kp1 
                        Ma_t(k,k) = -1; Ma_t(k,kp1) = 1;
                end
                if ii ~= k && ii ~= kp1 && ii <= M
                    if type_under(i) == -1 % type II
                        Ma_t(k,k) = -s; Ma_t(k,kp1) = 1;  
                        %if kp1 ~= ii
                         
                            Ma_t(k,ii) = s-1;
                        %end
                    else
                        Ma_t(k,k) = 1; Ma_t(k,kp1) = -s;  
                        %if kp1 ~= ii
                        
                            Ma_t(k,ii) = s-1;
                        %end
                    end
                end
                if ii > M
                    if type_under(i) == -1 % type II
                        Ma_t(k,k) = -t; Ma_t(k,kp1) = 1;  
                        %if kp1 ~= ii
                            Ma_t(k,ii) = s-1;
                        %end
                    else
                        Ma_t(k,k) = 1; Ma_t(k,kp1) = -t;  
                        %if kp1 ~= ii
                            Ma_t(k,ii) = s-1;
                        %end
                    end
                end
            end

            if k == M && M == 1 && ii > M
                    Ma_t(k,k) = 1-t; Ma_t(k,ii) = s-1;
            end

            if k > M && N>M+1
                if ii ==  k || ii == kp1
                    Ma_t(k,k) = -1; Ma_t(k,kp1) = 1;
                end

                if ii ~= k && ii ~= kp1 && ii > M
                    if type_under(i) == -1 % type II
                        Ma_t(k,k) = -t; Ma_t(k,kp1) = 1;  
                        %if kp1 ~= ii
                            Ma_t(k,ii) = t-1;
                        %end
                    else
                        Ma_t(k,k) = 1; Ma_t(k,kp1) = -t;  
                        %if kp1 ~= ii
                            Ma_t(k,ii) = t-1;
                        %end
                    end
                end

                if ii <= M 
                    if type_under(i) == -1 % type II
                        Ma_t(k,k) = -s; Ma_t(k,kp1) = 1;  
                        %if kp1 ~= ii
                            Ma_t(k,ii) = t-1;
                        %end
                    else
                        Ma_t(k,k) = 1; Ma_t(k,kp1) = -s;  
                        %if kp1 ~= ii
                            Ma_t(k,ii) = t-1;
                        %end
                    end
                end
            end

            if k == N && N == M+1 && ii <= M
                    Ma_t(k,k) = 1-s; Ma_t(k,ii) = t-1;
            end
        end

        % Knot invariant 
        Ma_t11 = Ma_t(2:N, 2:N);
        if M >= 1
            detAlex = abs((det(Ma_t11)/(s-1)));
            %f = factor(detAlex);
            
        end
        
    end
    
    
else
    detAlex = 0;
end

end

function y=pivot_single(v,i1,di,theta)
n=di;
i2=size(v,2);
J=[0 -n(3) n(2); n(3) 0 -n(1); -n(2) n(1) 0];
R=eye(3,3)+J*sin(theta)+J*J*(1-cos(theta));

v0=v(:,i1);
v1=R*(v(:,i1+1)-v0)+v0;
u1= v(:,i1+1);
y= [v1 (v(:,i1+2:i2)-u1)+v1];
end

function E=Ebend_vector(x,k_b)
d=diff(x,1,2);
thetas = acos(sum(d(:,1:end-1).*d(:,2:end)));
E= k_b * (thetas*thetas');
end
function y=CR(v,i1,i2,theta) % Crankshaft rotation
n=(v(:,i2)-v(:,i1))/norm(v(:,i2)-v(:,i1),2);
J=[0 -n(3) n(2); n(3) 0 -n(1); -n(2) n(1) 0];
R=eye(3,3)+J*sin(theta)+J*J*(1-cos(theta));
v0=repmat(v(:,i1),[1 i2-i1-1]);
y=R*(v(:,i1+1:i2-1)-v0)+v0;
end
function v = reptation(v1,i1,i2,c)
    v = v1;
    if c >= 0.5  % reptation backward
        v(:, i1:i2-1) = [v1(:,i1+1:i2-1)+repmat(v1(:,i1-1)-v1(:,i1),1,i2-i1-1) v1(:,i1-1)+v1(:,i2)-v1(:,i1)];
    else % reptation forward
        v(:, i1+1:i2) = [v1(:,i2+1)-v1(:,i2)+v1(:,i1) v1(:,i1+1:i2-1)+repmat(v1(:,i2+1)-v1(:,i2),1,i2-i1-1)];
    end
end
function energy = E_coulomb4_vector(r1, r2, d0,seg)
r = [r1, r2];
energy =0;
radius = d0/seg;
r_center = 0.5*(r(:,1:end-1) + r(:,2:end)); % 2x(N-1) 
r_center(:, length(r1)) = [];
ptsc =r';
D = pdist(ptsc);
testarr = sum(D <=radius);
    if testarr >0
       energy= realmax; % hard-wall
    end
    if energy ==0
        ptsc =r_center';
        D = pdist(ptsc);
        testarr = sum(D <=radius);

        if testarr >0
           energy= realmax; % hard-wall         
        end
    end
end
function energy = E_debye_huckel_no_neighbor_cubic(r1, r2, seg, rcutoff, fdinterpolation,choosearr,index2)

%r = [r1, r2];
r = ([r1,fliplr(r2)]);

% r(:, length(r1))= r1;
r_center = 0.5*(r(:,1:end-1) + r(:,2:end)); % 2x(N-1)     
%r_center(:, length(r1)) = [];

ptsc =r_center';
Darr = pdist(ptsc);

iarr = 2:size(r_center,2);
jarr = 1:size(r_center,2)-1;
test1 = (iarr-1).*((size(r,2))-iarr/2)+jarr-iarr+2;
test1 = [1,test1(1:end-1)];
test0 = 0:size(r_center,2)-2;
test1 = test1 - test0;
Darr(test1) =[];
Darr(index2) =[];

index = (Darr <= rcutoff/seg);
rhodist = seg * Darr(index);
%choosearr = nchoosek(1:length(r)-1,2);
%choosearr(test1,:) =[];

p1 = r(:,choosearr(index,1)); 
p2 = r(:,choosearr(index,1)+1); 
p3 = r(:,choosearr(index,2));
p4 = r(:,choosearr(index,2)+1); 
ei = p2-p1;
ej = p4-p3;
rij = seg*((p4+p3)/2 - (p2+p1)/2);
g1=dot(ei,rij)./rhodist;
g2=-dot(ej,rij)./rhodist;
%sig_num = dot(cross(ei,rij),cross(ej,rij));
sig = dot(cross(ei,rij),cross(ej,rij))./(vecnorm(cross(ei,rij)).*vecnorm(cross(ej,rij)));
sig(isnan(sig)) = 0;
% if (~isreal([acos(g1), acos(g2), acos(sig)]))
%     flag= 0;
%     allgvec = [g1,g2,sig];
%     allvec = [acos(g1), acos(g2), acos(sig)];
%     avec = abs(imag( allvec ));
%     [M,Ind] = max(avec);
%     disp( M )
%     disp( allvec(Ind) )
%    disp( allgvec(Ind) )
    %energy = (4.0453*sum(fdinterpolation(rhodist,real(acos(g1)),real(acos(g2)),real(acos(sig)))));
energy = (4.0453*sum(fdinterpolation(rhodist,real(acos(g1)),real(acos(g2)),real(acos(sig)))));
%else
%    energy = (4.0453*sum(fdinterpolation(rhodist,acos(g1),acos(g2),acos(sig))));
%end

end


