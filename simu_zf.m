clear ; 
% after default, the price of debt is recovery rate ->rr(st,bpol). In each period
% after default, the economy has prob mu to renegotiate (get bpr) and get reentry. 

load fiscal_zf2.mat

format compact

Sim = 500;
nos = 3000;
T = 3500;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% vectors in t = 1:T loop 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Z = zeros(T,1);    Y = zeros(T,1);    C = zeros(T,1);  
H = zeros(T,1);    Bad = zeros(T,1);    Q = zeros(T,1);    
B = zeros(T,1);    D = zeros(T,1);    % default with good record
PM = zeros(T,1);   Gov = zeros(T,1);
FE = zeros(T,1);   F = zeros(T,1) ;   % transfers payment
Tax = zeros(T,1);  % labour tax income
Re = zeros(T,1);   % indicator for re-entry
Cut = zeros(T,1); % haircut ratio
r_annual = ((1+rbase)^4-1)*100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% vectors in s = 1:sim loop 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ZZ = zeros(nos,Sim);    YY = zeros(nos,Sim);    CC = zeros(nos,Sim);  
BadB = zeros(nos,Sim);  DD = zeros(nos,Sim);    PMPM = zeros(nos,Sim);   
BB = zeros(nos,Sim);    DY = zeros(nos,Sim);    TBY = zeros(nos,Sim);    
HH = zeros(nos,Sim);    GG = zeros(nos,Sim);    FEFE = zeros(nos,Sim);  
FF = zeros(nos,Sim);
TaxT = zeros(nos,Sim);    CUT = zeros(nos,Sim);  

%
stats1 = zeros(Sim,5);    stats2 = zeros(Sim,8);  
stats3 = zeros(Sim,8);    stats4 = zeros(Sim,6);
st_govt = zeros(Sim,6);   st_tax = zeros(Sim,7);
st_transfer = zeros(Sim,7); 
st_cut = zeros(Sim,4);    st_B = zeros(Sim,6);  

% codegen simu_func -args {pdfz,T} -o simu_func_mex2
% zpath = simu_func_mex(pdfz,T) ; % generate shocks
% fpath = simu_func_mex(pdff,T) ;


for s = 1:Sim
    
  bad = 0;
  ib = floor(nb/3);
  reentry = rand(T,1);
  zpath = simu_func_mex2(pdfz,T)  ;
  fpath = simu_func_mex2(pdff,T) ; 
  
  for t = 1:T
    
    iz = zpath(t);
    ife = fpath(t);
    is = (ife-1)*nz + iz ;
    
    FE(t) = f(ife);
    B(t) = b(ib);
    ir = reentry(t);  %for re-entry

    if (bad==0) && (default(is,ib)==0) % 1. good record, choose to repay debt
        Z(t) = z(iz);
        D(t) = 0;
        Bad(t) = 0;
        pol = bp(is,ib);
        Q(t) = q(is,pol);
        tr = psi0_t + psi1_t*exp(z(iz)) ;
        c0 = 4/theta/(1-psi1_g+psi1_f)*(1-tr)*exp(z(iz))^2 ;
        ss = - q(is,pol)*(b(pol)-(1-eta)*b(ib)) + (eta + (1-eta)*coup+psi2_g)*b(ib)...
            + (psi0_g - psi0_f) - f(ife) ;
        C(t) = ( - ss + sqrt(ss^2 + c0))/2;
        H(t) = ( (1-tr)/theta/(1-psi1_g+psi1_f) * (1 + ss/C(t)) )^(1/(chi+1)) ;
        Y(t) = exp(z(iz))*H(t) ;
        Tax(t) = tr*Y(t);
        F(t) = psi0_f + psi1_f*Y(t) + f(ife) ;
        Gov(t) = psi0_g + psi1_g*Y(t) + psi2_g*b(ib) ;
        Re(t) = 0;
    end

    if (bad==0) && (default(is,ib)==1) % 2. good record, choose to default (pol=id) and fall into autarky
        Z(t) = za(iz);
        D(t) = 1;
        Bad(t) = 0;
        pol = ib;  % no new borrowing and no haircut/bargaining
        Q(t) = rr(is,pol);
        tr = psi0_t + psi1_t*exp(za(iz)) ;
        ca0 = 4/theta/(1-psi1_g+psi1_f) * (1-tr) * exp(za(iz))^2 ;
        sa = (psi0_g - psi0_f) - f(ife) ;
        C(t) = ( - sa + sqrt(sa.^2 + ca0))/2;
        H(t) = ( (1-tr)/theta/(1-psi1_g+psi1_f) * (1 + sa/C(t)) )^(1/(chi+1)) ;
        Y(t) = exp(za(iz))*H(t) ;
        Tax(t) = tr*Y(t);
        F(t) = psi0_f + psi1_f*Y(t) + f(ife) ;
        Gov(t) = psi0_g + psi1_g*Y(t) + psi2_g*b(ib) ;
        Re(t) = 0; 
    end

    if (bad==1) && (ir>mu) % 3. bad record, autarky and no bargaining
        Z(t) = za(iz);
        D(t) = 0;
        Bad(t) = 1;
        pol = ib;  % no new borrowing and no haircut/bargaining
        Q(t) = rr(is,pol);
        tr = psi0_t + psi1_t*exp(za(iz)) ;
        ca0 = 4/theta/(1-psi1_g+psi1_f)* (1-tr)*exp(za(iz))^2 ;
        sa = (psi0_g - psi0_f) - f(ife) ;
        C(t) = ( - sa + sqrt(sa.^2 + ca0))/2;
        H(t) = ( (1-tr)/theta/(1-psi1_g+psi1_f) * (1 + sa/C(t)) )^(1/(chi+1)) ;
        Y(t) = exp(za(iz))*H(t) ;
        Tax(t) = tr*Y(t);
        F(t) = psi0_f + psi1_f*Y(t) + f(ife) ;
        Gov(t) = psi0_g + psi1_g*Y(t) + psi2_g*b(ib) ;
        Re(t) = 0;
    end

    if (bad==1) && (ir<=mu) && (default(is,ib)==0) % 4. bad record, reenter, bargain and not to default
        Z(t) = z(iz);
        D(t) = 0;
        Bad(t) = 0;
        pol = bpr(is);   % new debt outstanding, after haircut
        Q(t) = q(is,pol);
        tr = psi0_t + psi1_t*exp(z(iz)) ;
        c0 = 4/theta/(1-psi1_g+psi1_f)*(1-tr)*exp(z(iz))^2 ;
        ss = - q(is,pol)*(b(pol)-(1-eta)*b(ib)) + (eta + (1-eta)*coup+psi2_g)*b(ib)...
            + (psi0_g - psi0_f) - f(ife) ;
        C(t) = ( - ss + sqrt(ss^2 + c0))/2;
        H(t) = ( (1-tr)/theta/(1-psi1_g+psi1_f) * (1 + ss/C(t)) )^(1/(chi+1)) ;
        Y(t) = exp(z(iz))*H(t) ;
        Tax(t) = tr*Y(t);
        F(t) = psi0_f + psi1_f*Y(t) + f(ife) ;
        Gov(t) = psi0_g + psi1_g*Y(t) + psi2_g*b(ib) ;
        Re(t) = 1;  
        Cut(t) = 1 - b(bpr(is))/b(bp(is,ib));
    end

    if (bad==1) && (ir<=mu) && (default(is,ib)==1) % 5. bad record, reenter, bargain and to default (not happy with bargaining)
        Z(t) = za(iz);
        D(t) = 0;
        Bad(t) = 1;
        pol = ib;  % new debt outstanding, after haircut
        Q(t) = rr(is,pol);
        tr = psi0_t + psi1_t*exp(za(iz)) ;
        ca0 = 4/theta/(1-psi1_g+psi1_f) * (1-tr)* exp(za(iz))^2 ;
        sa = (psi0_g - psi0_f) - f(ife) ;
        C(t) = ( - sa + sqrt(sa.^2 + ca0)) / 2 ;
        H(t) = ( (1-tr)/theta/(1-psi1_g+psi1_f) * (1 + sa/C(t)) )^(1/(chi+1)) ;
        Y(t) = exp(za(iz))*H(t) ;
        Tax(t) = tr*Y(t);
        F(t) = psi0_f + psi1_f*Y(t) + f(ife) ;
        Gov(t) = psi0_g + psi1_g*Y(t) + psi2_g*q(is,pol) ;
        Re(t) = 0;  
    end
    
    ib = pol;
    PM(t) = (((eta+(1-eta)*coup)/Q(t)+1-eta )^4-1)*100 - r_annual; % country risk premium
    bad = Bad(t) + D(t); 
    if bad > 1; error("Credit record wrong"); end   
  end

dy1 = B./Y*100;       
nx = (B(2:T)-(1-eta)*B(1:T-1)).*Q(1:T-1) - (eta+(1-eta)*coup)*B(1:T-1);

%%%%%% eliminate burning period %%%%%% 
ZZ(:,s) = exp(Z(T-nos:T-1));     
TBY(:,s) = nx(T-nos:T-1)./Y(T-nos:T-1);
CC(:,s) = C(T-nos:T-1);   
HH(:,s) = H(T-nos:T-1); 
BB(:,s) = B(T-nos:T-1);            
DY(:,s) = dy1(T-nos:T-1);  
YY(:,s) = Y(T-nos:T-1);

GG(:,s) = Gov(T-nos:T-1);
FF(:,s) = F(T-nos:T-1);
FEFE(:,s) = FE(T-nos:T-1);
TaxT(:,s) = Tax(T-nos:T-1);
BadB(:,s) = Bad(T-nos:T-1);
DD(:,s) = D(T-nos:T-1);   
PMPM(:,s) = PM(T-nos:T-1); 

x = find(BadB(:,s)==0 & DD(:,s)==0);
def_frequency = mean(DD(:,s))*4*100; %default frequency per century
zz = ZZ(:,s);     zz = zz(x);           
yy = YY(:,s);     yy = yy(x);
cc = CC(:,s);     cc = cc(x);           
hh = HH(:,s);     hh = hh(x);
dd = BB(:,s);     dd = dd(x);                     
dy = DY(:,s);     dy = dy(x);  
tby = TBY(:,s);   tby = tby(x);         
pm = PMPM(:,s);   pm = pm(x);
gov = GG(:,s);    gov = gov(x);
ff = FF(:,s) ;    ff = ff(x) ;
fe = FEFE(:,s);   fe = fe(x);
tax = TaxT(:,s);  tax = tax(x);
% yy_hp = YY_hp(x);

[~,yy_hp] = hpfilter(log(yy),1600);
[~,cc_hp] = hpfilter(log(cc),1600);
[~,gov_hp] = hpfilter(log(gov),1600);
[~,tax_hp] = hpfilter(log(tax),1600);
[~,ff_hp] = hpfilter(log(ff),1600);

kk = Re==1;  haircut = Cut(kk)*100;

stats1(s,:)= [def_frequency mean(dy) mean(pm) std(pm) std(yy_hp) ] ;

stats2(s,:)= [std(cc_hp)/std(yy_hp) std(tby)/std(yy_hp) std(zz) corr(zz(1:end-1),zz(2:end))...
    corr(cc,zz) corr(pm,zz) corr(pm,dy) corr(dy,yy)] ;  

stats3(s,:)= [corr(yy,cc) corr(yy,tby) corr(yy(1:end-1),yy(1:end-1)) corr(cc(1:end-1),cc(2:end)) ...
    corr(yy,zz) corr(pm,cc) corr(pm,yy) corr(pm,tby)] ;

stats4(s,:)= [corr(hh,yy), corr(hh(1:end-1),hh(2:end)), corr(hh,pm), ...
    corr(zz,hh), std(hh), corr(tby(1:end-1),tby(2:end)) ];

st_govt(s,:)= [mean(gov./yy), std(gov./yy), corr(gov_hp,yy_hp),  ...
    corr(gov./yy,dy), std(gov_hp), corr(gov(1:end-1),gov(2:end))] ;

st_tax(s,:)= [mean(tax./yy), std(tax./yy), corr(tax_hp,yy_hp), ...
    corr(tax(1:end-1),tax(2:end)), corr(tax,pm), corr(tax./yy,pm), ...
    corr(tax./yy,dy) ] ;

st_transfer(s,:) = [mean(ff./yy), std(ff./yy), corr(log(ff),log(yy)), ...
    corr(ff(1:end-1),ff(2:end)), corr(ff,pm), corr(ff./yy,pm), ...
    corr(ff./yy,dy) ] ;

st_cut(s,:)=[mean(haircut) median(haircut) min(haircut) max(haircut)];

st_B(s,:)=[max(B), median(B), min(B), max(dd), median(dd), min(dd)];

end

%%
moment1 = median(stats1)
moment2 = median(stats2)
moment3 = median(stats3)
moment4 = median(stats4)
moment_govt = median(st_govt)
moment_tax = median(st_tax)
moment_transfer = median(st_transfer)
moment_haircut = median(st_cut)
max_Bgrid = median(st_B)
mean_Y = median(mean(YY));

%%
lab = zeros(Sim,nb);
for s = 1:Sim
    % x = find( BB(:,s)==0 & FF(:,s)==0 & DD(:,s) < b(end));
    x = find( BadB(:,s)==0 & DD(:,s)==0 & BB(:,s) );
    d1 = BB(:,s);     
    d1 = d1(x);
    for i = 1:nb
    lab(s,i) = mean(d1 == b(i));
    end
end

lab_mean = mean(lab);

%%
figure
plot(b,lab_mean,'x-' ,'linewidth', 2)
xlabel('b')
ylabel('lab')

%%
save simu_zf2.mat moment1 moment2 moment3 moment4 moment_govt moment_tax moment_transfer...
    moment_haircut max_Bgrid lab_mean b mean_Y Sim moment_tax

