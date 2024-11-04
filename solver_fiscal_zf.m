function [bp,bpr,default,probdef,q,rr,vp,vd,za] = solver_fiscal_zf(b,f,z,pdff,pdfz,para)

mu = 0.0385;  % prob of re-entry 
chi = 1;  % inverse of Frisch elasticity
rbase = 0.01; % baseline interest rate
coup = 0.0125; % long-term bond, coupon rate
eta = 0.035; % 1/eta quarters, long-term bond, average maturity

%% release calibrated parameters 
% para = [phi0; phi1; alfa; betta; psi0_f; psi1_f; psi0_g; psi1_g; psi2_g; psi0_t;
%     psi1_t; theta; sigg_bp; sigg_defp; sigg_bpr ];

phi0 = para(1);      phi1 = para(2);     alfa = para(3);       betta = para(4); 
psi0_f = para(5);    psi1_f = para(6);
psi0_g = para(7);    psi1_g = para(8);   psi2_g = para(9) ;   
psi0_t = para(10);   psi1_t = para(11); 
theta = para(12);    sigg_bp=para(13);   sigg_defp=para(14);
sigg_bpr=para(15);  

%%
nz = length(z);
nf = length(f);
nb = length(b);
ns = nz*nf; % number of grids for exogenous shocks

za = exp(z) - max(0, -phi0*exp(z) + phi1*exp(z).^2);  % TFP loss function
za = log(za);

pdf = kron(pdff,pdfz); %Joint distribution, each g includes all t grids, is = (ig-1)*nt+it
pdf = sparse(pdf./(((sum(pdf,2)))*ones(1,ns)));  % norm(pdf-pdf0)

% pdf = pdf./(((sum(pdf,2)))*ones(1,ns));
%  ig=2; it=2; iz =1; is = ((ig-1)*nt+it-1)*nz+iz
%  pdf1 = kron(pdfg(ig,:),pdft(it,:)); pdf1 = kron(pdf1,pdfz(iz,:)); norm(pdf(is,:) - pdf1)

% ge1 = zeros(ns,1); zaza1 = ge1; tt1 = ge1; zz1 =ge1;
% for ig = 1:ng
%   for it = 1:nt
%     for iz = 1:nz
%       is = ((ig-1)*nt+it-1)*nz+iz;
%       ge1(is) = g(ig);
%       tt1(is) = tau(it);
%       zaza1(is) = za(iz);
%       zz1(is) = z(iz);
%     end
%   end
% end

zz = repmat(z,nf,1);   % norm(zz - zz1)
zaza = repmat(za,nf,1);  % norm(zaza - zaza1)
fe = repmat(f',nz,1); fe = fe(:);  % norm(ge1 - ge)
tt = psi0_t + psi1_t*exp(zz) ;

%% consumption and labour
% to solve c* = (1/2)*(-st + sqrt(st^2 + 4/theta*(1-psi1_g+psi1_f)*(1-tauh)*exp(z)^2) ), 
% where s = - q[b' - (1-eta)b - psi2_g] + [eta + (1-eta)coup]b + (psi0_g - psi0_f) - fe

c0 = (4/theta*(1-psi1_g+psi1_f))*(1-tt).*exp(zz).^2 ; 
% c0 is the 4/theta/(1-psi1_g+psi1_f)*(1+tauc)*(1-tauh)*exp(z)^2 component of c*

%%%%% for autarky, caut = (1/2)*(-saut + sqrt(saut^2 + 4/theta*(1-psi1_g+psi1_f)*(1+tauc)*(1-tauh)*exp(za)^2) ), 
% where saut = (psi0_g - psi0_f) - fe
ta = psi0_t + psi1_t*exp(zaza) ;
saut = (psi0_g - psi0_f) - fe ;
caut = 4/theta*(1-psi1_g+psi1_f)*(1-ta).*exp(zaza).^2 ; 
caut = (- saut + sqrt(saut.^2 + caut)) / 2 ;
% h = ( (1 + s/c)/(1-psi1_g+psi1_f)/theta ) ^(1/(1+chi))  
ua = log(caut) - (1 + saut./caut).*(1-ta)/(chi+1)/(1-psi1_g+psi1_f) ;

dist = 1;    vaut = zeros(ns,1); %value of autarky
while dist > 1e-8
    vautnew = ua + betta*pdf*vaut; 
    dist = max(abs(vautnew(:)-vaut(:)));
    vaut = vautnew;
end  

% save fiscal_gzt_vaua.mat ua vaut za
% clear vautnew ca sa pdfg pdft pdfz zaza zz te ta

%% to incorporate taste shocks
epsi = 10e-16;
cv_bpr = sigg_bpr*log(epsi);  % critical value
probDcre = zeros(ns,1);
probVp = zeros(ns,1);
cv_bp = sigg_bp*log(epsi); % critical value

%% Initialize the Value functions
V = zeros(ns,nb);  %continue repaying
vp = zeros(ns,nb); %the value of good standing
vpnew = vp;
vd = zeros(ns,1); %value of default
vdnew = vd;
ev = vp;

bp = zeros(ns,nb); %debt policy function (expressed in indices)  
bpr = zeros(ns,1); % debt policy (index) when decided renegotiate (right after every default)
q = ones(ns,nb)/(1+rbase); %initial price of debt
qnew = zeros(ns,nb);
rr = 0.5*ones(ns,nb)/(1+rbase);
probdef = zeros(ns,nb);

%% Execute the VFI
diff = 1;    its = 1;
tol = 1e-6;    maxits = 1400;
smctime = tic;  totaltime = 0;

while diff > tol && its< maxits

ev = betta*pdf*V;

  parfor is = 1:ns

      W = zeros(1,nb,'double') ;
      probbp = zeros(nb,nb);
      qv = q(is,:) ;
      tax = tt(is) ;

      for ib = 1:nb

          bib = b(ib) ;

          for ibp = 1:nb

              s = (eta+(1-eta)*coup+psi2_g) * bib - (b(ibp)-(1-eta)*bib).*qv(ibp) ...
              + (psi0_g - psi0_f) - fe(is);
              c = ( -s + sqrt(s.^2 + c0(is)) ) / 2 ;
              W(ibp) = log(c) - (1-tax)*(1 + s./c)/(1+chi)/(1-psi1_g+psi1_f) + ev(is,ibp);

          end
          vpnew(is,ib) = max( W ) ;
          indix = W - vpnew(is,ib) - cv_bp > 0 ; 
          theExp = exp( (W(indix)- vpnew(is,ib)) / sigg_bp ) ;
          probbp(ib,indix) = theExp ./ sum( theExp ) ; 
      end

    qnew(is,:) =  (eta + (1-eta)*(coup + sum(probbp.*qv,2) )) ;

  end

  parfor is = 1:ns
    
    vpnewv = vpnew(is,:) ;
    %%%%% bargaining of debt restructure
    Dsov = max( 0, vpnew(is,:) - vaut(is) );
    Dcre = (eta + (1-eta)*(coup + q(is,:))).*b';
    W_is = Dsov.^alfa.*Dcre.^(1-alfa);

    bargain_is = max(W_is,[],2);

    indix = W_is - bargain_is - cv_bpr > 0;
    probbpr = 1 ./ sum( exp( ( W_is(indix)'-W_is(indix))/sigg_bpr ) ); 
    probDcre(is) = sum(probbpr.*Dcre(indix));
    probVp(is) = sum(probbpr.*vpnewv(indix));

  end
  
  vdnew = ua + betta*pdf*( mu*probVp + (1-mu)*vd );
  probdef = 1./( 1 + exp((vpnew - vdnew)/sigg_defp) ); % prob of default

  rrnew = pdf*((1-mu)*rr + mu*probDcre./b')/(1 + rbase); 
%rr: recovery rate(Yue,2010), the value of defaulted b is reduced to chi percent of the unpaid b
  qnew = pdf*(probdef.*rrnew + (1-probdef).*qnew)/(1 + rbase);
  
  diff = max(abs(qnew(:)-q(:))) + max(abs(rrnew(:)-rr(:))) + max(abs(vpnew(:)-vp(:))) ;
  % diff = max(abs(qnew(:)-q(:))) + max(abs(vpnew(:)-vp(:))) ;

  q = qnew;
  vp = vpnew;
  vd = vdnew;
  rr = rrnew;
  V = max(vp,repmat(vd,1,nb));

totaltime = totaltime + toc(smctime);
avgtime   = totaltime/its;

  if mod(its, 50) == 0
    fprintf('%4.0f ~%4.7f ~%4.5fs ~%4.5fs \n', its, diff, totaltime, avgtime);
  end

its = its+1;
smctime = tic; % re-start clock

end

default = vpnew < repmat(vdnew,1,nb);


%%%% find bp and bpr
parfor is = 1:ns

      W = zeros(1,nb,'double') ;
      tax = tt(is) ;

      for ib = 1:nb

          bib = b(ib) ;

          for ibp = 1:nb
              s = (eta+(1-eta)*coup+ psi2_g) * bib - ( b(ibp) - (1-eta)*bib).*q(is,ibp) ...
              + (psi0_g - psi0_f) - fe(is);
              c = ( -s + sqrt(s.^2 + c0(is)) ) / 2 ;
              W(ibp) = log(c) - (1-tax)*(1 + s./c)/(1+chi)/(1-psi1_g+psi1_f) + ev(is,ibp);
          end
          [~, bp(is,ib)] = max(W,[],2);
      end
end

parfor is = 1:ns

    %%%%% bargaining of debt restructure
    Dsov = max( 0, vpnew(is,:) - vaut(is) );
    Dcre = (eta + (1-eta)*(coup + q(is,:))).*b';
    W_is = Dsov.^alfa.*Dcre.^(1-alfa);

    [~, bpr(is,:)] = max(W_is,[],2);

end

end


