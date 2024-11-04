clear ; 
% after default, the price of debt is recovery rate ->rr(st,bpol). In each period
% after default, the economy has prob mu to renegotiate (get bpr) and get reentry. 

load fiscal_zf2.mat

format compact

rng(220)

nos = 1000000;
T = 1010000;

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
Z = exp(Z(T-nos:T-1));     
TBY = nx(T-nos:T-1)./Y(T-nos:T-1);
C = C(T-nos:T-1);   
H = H(T-nos:T-1); 
B = B(T-nos:T-1);            
BY = dy1(T-nos:T-1);  
Y = Y(T-nos:T-1);
Gov = Gov(T-nos:T-1);
F = F(T-nos:T-1);
FE = FE(T-nos:T-1);
Tax = Tax(T-nos:T-1);
Bad = Bad(T-nos:T-1);
D = D(T-nos:T-1);   
PM = PM(T-nos:T-1); 

% [~,yy_hp] = hpfilter(log(Y),1600);
% [~,cc_hp] = hpfilter(log(C),1600);
% [~,gov_hp] = hpfilter(log(Gov),1600);
% [~,tax_hp] = hpfilter(log(Tax),1600);
% [~,ff_hp] = hpfilter(log(F),1600);

%%
lab = zeros(1,nb);
x = find( Bad ==0 & D ==0 & B );
d1 = B;     
d1 = d1(x);
for i = 1:nb
    lab(i) = mean(d1 == b(i));
end
figure
plot(b,lab,'x-' ,'linewidth', 2)
xlabel('b')
ylabel('lab')

%%
% save simu_longpath.mat Y C Gov Tax F FE PM FE Bad Z TBY D
