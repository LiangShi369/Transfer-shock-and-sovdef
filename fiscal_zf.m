clear;

format compact
%% parameters to calibrate
phi0 = 0.3702 ; % higher abs(phi0) => the initial first grid of zaut<z is at higher zi 
             %  (smaller z index) so higher default
             % if dy too high, increase abs(phi)
phi1 = 0.403 ; % higher phi1 => yaut is further lower than y, so fewer default
              % and increases E(d/y)
alfa = 0.66 ;  % borrower's bargaining power, higher alfa => higher haircut & higher default
betta = 0.977 ; % discount factor, to match mean(rs)

%% Gov transfers rule
% f_t = psi0 + psi1*y_t
psi0_f = -0.115 ; % to match mean(f/y) = 14.7% 
psi1_f = 0.26 ; % corr(ff,yy) = 0.26

%% Gov spending rule:
% g_t = psi0 + pis1*y_t + psi2*b_t
psi0_g = -0.07 ; % to match mean(g/y)
psi1_g = 0.27 ;  % from an estimation of fiscal rule, in the Data folder
psi2_g = 0.004 ; 

%% tax rule
psi0_t = 0.128;     % to match mean(t/y)
psi1_t = -0.0526; % from an estimated fiscal rule, in the Data folder

%% parameters for state variables
nz = 35;   % # grids of TFP                  21
nf = 35;   % # grids of govt spending shock  21
nb = 300;  % # of grid points for debt   150

rhoz = 0.90;   
sdz = 0.0152 ;  % stdev of log TFP
rhof = 0.95;    % autocorrelation of g^e process
sdf = 0.007;  % error term of g^e  
bupper = 1 ;   % debt grid 

std_z = sdz/sqrt(1-rhoz^2) ;
std_f = sdf/sqrt(1-rhof^2) ;

%% parameters for taste shocks
sigg_bpr = 0.0003;  % for renegotiation
sigg_bp = 0.0003;   % for debt policy in good status
sigg_defp = 0.0003; % for default decision

%% pre-set parameters
mu = 0.04;  % prob of re-entry 
chi = 1;  % inverse of Frisch elasticity
rbase = 0.01; % baseline interest rate
coup = 0.0125; % long-term bond, coupon rate
eta = 0.035; % 1/eta quarters, long-term bond, average maturity

% endowment grids
width = 4.5;  % 3.7 (two-tail) sigma both sides covers 99.98% of the distribution
[z,pdfz] = tauchen(nz,0,rhoz,sdz,width);
[f,pdff] = tauchen(nf,0,rhof,sdf,width);

theta = (1 - psi0_t - psi1_t) /( 1 - psi0_g - psi1_g + psi0_f +psi1_f ) ;

% debt grids
blower = 0.3;
b = blower:(bupper-blower)/(nb-1):bupper; 
b = b(:); 

para = [phi0; phi1; alfa; betta; psi0_f; psi1_f; psi0_g; psi1_g; psi2_g; psi0_t;
    psi1_t; theta; sigg_bp; sigg_defp; sigg_bpr ];

% [bp,bpr,default,probdef,q,rr,vgood,vbad,za] = solver_fiscal_zf(b,f,z,pdff,pdfz,para);
%%
cfg = coder.config('mex');
cfg.GenerateReport = true;

tic
codegen -config cfg solver_fiscal_zf -args {b,f,z,pdff,pdfz,para} -o solver_fiscal_zf_mex
toc

[bp,bpr,default,probdef,q,rr,vgood,vbad,za] = solver_fiscal_zf_mex(b,f,z,pdff,pdfz,para);

%%
vdiff = vgood - vbad;
vbad = repmat(vbad,1,nb);

zz = repmat(z,nf,1);
fe = repmat(f',nz,1); fe = fe(:);

save fiscal_zf2.mat alfa b betta bp bpr chi coup default eta f psi0_g psi1_g psi2_g...
    mu psi0_t psi1_t psi0_f psi1_f nf nz nb pdff pdfz phi0 phi1 q rbase rr ...
    rhof rhoz sigg_bp sigg_bpr sigg_defp sdf sdz theta vbad vgood width z za ...
    probdef

%% analyse

% iz = 13 ;
% ib = 220 ;
% 
% ife_low = 4 ; 
% is1 = (ife_low-1)*nz + iz ;
% 
% ife_high = 32 ;
% is2 = (ife_high-1)*nz + iz ;
% 
% [q(is1,bp(is1,ib)) q(is2,bp(is2,ib)); default(is1,bp(is1,ib)) default(is2,bp(is2,ib))]  % [low f, high f ]

%%

