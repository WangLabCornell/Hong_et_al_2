function [vj1,vj2,dLk]= braiding_init_configs(n0, dtopbp, dbotbp, dLk0,seg, r)

toln = 1e-10;
rise = 0.338;       % nm
L0 = n0*rise;       % contour length in nm for both strands
n = round(L0/seg);  % number of segments, for both strand, n is even
if mod(n,2) ~= 0
    n = n+1;      % if n is odd increase length by 1 segment of 10 nm
end

Zfar = 200;        % distance to infinity, in segment

%% distance at bottom, 2 chains are anchored along x, distance by 2*e
%% bottom anchoring of chain 2 is on +x compared to bottom anchoring of chain 1
et = (dtopbp*0.338)/2; % 'e' in nm
et = et/seg; % e in seg units

eb = (dbotbp*0.338)/2; % 'e' in nm
eb = eb/seg; % e in seg units

%% spring constant for additional energy terms

%for r= 3       % r is the radius of the braid, in unit of segment
etcheck =realmax;
ebcheck = realmax;
n2t = max(5,ceil(r-5))-1;    
n2b = n2t;
      n1 = n-n2b-n2t;

% if  mod(n1,2)==1
%     n1=n1-1;
%     n2t=n2t+1;
% end
if dLk0==1
    fact = 4;
else
    fact = 2;
end

t_inc = et>eb;

while abs(etcheck - et)>1e-1 || abs(ebcheck - eb)>1e-1
    if t_inc
        n2t = n2t + 2;    
        if ebcheck<eb
        n2b = n2b + 2;
        end
    else
        n2b = n2b + 2;
        if etcheck<et
            n2t = n2t + 2;    
        end
    end

if dLk0 > 0
    %% constructing the first conformation

    n1 = n-n2b-n2t;
  
    % if  mod(n1,2)==1
    %     n1=n1-1;
    %     n2t=n2t+1;
    % end
    
    phii = linspace(0, 2 * pi * dLk0, n1+1 );
    phii(end) = 2 * pi * dLk0;
    p =  2 *pi;
    v1_braid = [-r * cos(phii);  -r * sin(phii) ; p/2/pi * phii];
    v2_braid = [r * cos(phii);  r * sin(phii) ; p/2/pi * phii];
    d1 = diff(v1_braid,1,2); 
    v1_braid = v1_braid ./ mean(sqrt(sum(d1.^2 ,1))); % normalization
    v2_braid = v2_braid ./ mean(sqrt(sum(d1.^2 ,1))); % normalization

    alphab = real(asin((v1_braid(1,1)+eb)/n2b));
    alphat = real(asin((v1_braid(1,end)+et)/n2t));
    
    v_1_addlower = zeros(3,n2b+1);
    v_2_addlower = zeros(3,n2b+1);

    v_1_addupper = zeros(3,n2t+1);
    v_2_addupper = zeros(3,n2t+1);

    for j = 0 : n2b
        v_1_addlower(:,j+1) = [-eb + j * sin(alphab), 0, j * cos(alphab)]';
        v_2_addlower(:,j+1) = [eb - j * sin(alphab), 0, j * cos(alphab)]';

    end
    for j = 0 : n2t
        v_1_addupper(:,j+1) = [-et - j * sin(alphat), 0, j * cos(alphat)]';
        v_2_addupper(:,j+1) = [et + j * sin(alphat), 0, j * cos(alphat)]';
    end

    v1_braid_all = [v_1_addlower(:,1:end-1) v1_braid + repmat(v_1_addlower(:,end)-v1_braid(:,1), 1, size(v1_braid,2))];
    v2_braid_all = [v_2_addlower(:,1:end-1) v2_braid + repmat(v_2_addlower(:,end)-v2_braid(:,1), 1, size(v1_braid,2))];
    v1_braid_all = [v1_braid_all(:,1:end-1)  v_1_addupper + repmat(v1_braid_all(:,end)-v_1_addupper(:,1), 1, size(v_1_addupper,2))];
    v2_braid_all = [v2_braid_all(:,1:end-1)  v_2_addupper + repmat(v2_braid_all(:,end)-v_2_addupper(:,1), 1, size(v_2_addupper,2))];
else
    v = zeros(3,n+1);
    v(3,:) = 0:n;
    
    v1_braid_all = v; v1_braid_all = v1_braid_all + repmat([-eb,0,0]',1, n+1);
    v2_braid_all = v; v2_braid_all = v2_braid_all + repmat([et,0,0]',1, n+1);
end
ebcheck = abs((v1_braid_all(1,1) - v2_braid_all(1,1))/2);
etcheck = abs((v1_braid_all(1,end) - v2_braid_all(1,end))/2);
end
plot3(v1_braid_all(1,:),v1_braid_all(2,:),v1_braid_all(3,:), 'r')
hold on
plot3(v2_braid_all(1,:),v2_braid_all(2,:),v2_braid_all(3,:), 'b')
hold off
axis equal 

vj1= v1_braid_all;
vj2 = v2_braid_all;
dv_extend1 = [ -Zfar, vj1(2,end), vj1(3,end) ;-Zfar,0,0]';   % wings of a square
dv_extend2 = [  Zfar, vj2(2,end), vj2(3,end) ; Zfar,0,0]';   % wings of a square

v1_extend = [vj1 , dv_extend1 , vj1(:,1)];
v2_extend = [vj2 , dv_extend2 , vj2(:,1)];
v1_extendxz = [v1_extend(1,:); v1_extend(3,:); v1_extend(2,:)] ;
v2_extendxz = [v2_extend(1,:); v2_extend(3,:); v2_extend(2,:)] ;   
  
choices = nchoosek(1:length(v1_extendxz)+length(v2_extendxz)-1,2);


dLk1 = Alex_link_speedx1v3(v1_extendxz, v2_extendxz, -1 + toln, -1 + toln , choices);
dLk = Alex_link_speedx1v3(v1_extendxz, v2_extendxz, -1.1 + toln, -1.1 + toln,choices);

disp(['GIVEN : ' num2str(dLk0)])
disp(['GIVEN Radius: ' num2str(r)])
disp(['Del (t=-1.1) : ' num2str(dLk)])
disp(['Del (t=-1) : ' num2str(dLk1)])

disp(['Extension : ' num2str(mean([vj1(3,end) vj2(3,end)]))])
%disp(['Hard Wall Energy vector (d0=' num2str(d0) ') : ' num2str(E_coulomb4_vector(vj1,vj2,5,seg))])
%disp(['Repulsion Energy vector (d0=' num2str(d0) ') : ' num2str(E_debye_huckel_no_neighbor_cubic(vj1,vj2,seg, rcutoff, ftrial,choices_rep,index2))])

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