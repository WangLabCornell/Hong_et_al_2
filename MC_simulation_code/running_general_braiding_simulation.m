%%

% specify end distances a and b 
afact = 1;
bfact = 10;

twoe = 21; % in bp

two_et = afact * twoe; % in bp
two_eb = bfact* twoe; % in bp

% specify DNA length and its discretization
n0=2184; %bp
rise = 0.338;    
lc = n0*rise;  %nm
seg=5; %nm
n = round(lc/seg);

% specify mechanics - force and turns added 
f=2;% pN
dlkarr = round(linspace(1,27,11)); % 2 pN

% initialize geometry
config1 = zeros(3,n+1,length(dlkarr) );
config2 = zeros(3,n+1,length(dlkarr) );
dlk_input_arr = zeros(1, length(dlkarr) );


for si = 1:11
    [tempconfig1, tempconfig2, tempdLk] =  braiding_init_configs(n0, two_et, two_eb, dlkarr(si),seg,3); %get starting config chain_generating and slight wlc model
    config1(:,:,si) =  tempconfig1;
    config2(:,:,si) =  tempconfig2;
    dlk_input_arr(si) =tempdLk;
end

% load debyye-huckel interpolant, see Klenin,... Langowski BJ 1998 
load('interpolant_acos_separation_length0.64312nm_veff_5.0803epernm.mat')

%%
% specify MC parameters
d0=5;
nc = 1e6;
trials =100;
kmax=1e8;
repeat = 1;

% specify spring constants
Lp = 48.3; %nm
kT = 4.1;
kt = kT * 50;     % to align the direction of the two connecting vectors at bottom and at top
ko = kT * 200;    % to align the chain to +z
kd = kT * 20;     % to fix the distance between 2 top anchoring points

% start MC simulations
for si = 5
   ttconfig1 = config1(:,:,si);
   ttconfig2 = config2(:,:,si);
   [vj1, vj2] = braiding_simu_gen(n0, f, Lp, two_et, two_eb, dlkarr(si), d0, nc, trials, kmax, ftrial, ttconfig1, ttconfig2, seg, 225, kt, ko,kd,repeat)
end

%%
