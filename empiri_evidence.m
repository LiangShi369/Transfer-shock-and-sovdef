clear;  clc

format compact
% Greece data from 1999Q1 to 2018Q2

%% GDP current price

% gdpQ = readmatrix('data collection.xlsx','Sheet','y','Range','b27:b104'); 
% [~,y] = hpfilter(log(gdpQ),1600) ;
% 
% d_y = NaN(length(gdpQ),1);
% d_y(2:end) = gdpQ(2:end) - gdpQ(1:end-1);
% 
% transfer_exp = readmatrix('data collection.xlsx','Sheet','transfer','Range','c11:c88'); 
% transfer_rev = readmatrix('data collection.xlsx','Sheet','transfer','Range','o11:o88'); 
% transfer = transfer_exp./transfer_rev;  % transfer expense over transfers rev
% d_trans = NaN(length(transfer),1);
% d_trans(2:end) = transfer(2:end) - transfer(1:end-1);
% 
% 
% debt = readmatrix('data collection.xlsx','Sheet','gov debt','Range','c11:c88'); 
% [~,bb] = hpfilter(log(debt),1600) ;
% b = nan(length(debt),1);
% b(length(b) - length(bb)+1:end,1) = bb;
% 
% govc = readmatrix('data collection.xlsx','Sheet','gov c','Range','c27:c104'); 
% [~,gg] = hpfilter(log(govc),1600) ;
% g = nan(length(govc),1);
% g(length(g) - length(gg)+1:end,1) = gg;
% 
% spreads = readmatrix('data collection.xlsx','Sheet','rs','Range','b27:b104'); 
% spreads = spreads./100;
% d_rs = NaN(length(spreads),1);
% d_rs(2:end) = spreads(2:end) - spreads(1:end-1);
% 
% yields = readmatrix('data collection.xlsx','Sheet','r','Range','b11:b88');
% q = (eta + (1-eta)*coup)./( (1+yields/100).^0.25 + eta -1 ) ;
% 
% rs = hpfilter(spreads,1600) ;
% 
% tax = readmatrix('data collection.xlsx','Sheet','taxes amount','Range','h11:h88');
% 
% 
% 
% save data_transfers.mat spreads q govc debt gdpQ transfer transfer_exp transfer_rev tax

load data_transfers.mat

dum_eap = zeros(length(transfer),1) ;
dum_eap(47:end) = 1;

% default period 2012Q2(the 54th obvervation) to 2016Q4 (the 72th)
dum_predef = zeros(length(transfer),1) ;
dum_predef(1:53) = 1;

%% How does T_t positively impact spreads (during predef periods)? 

tb_rs_predef = table(log(transfer(1:53)), log(transfer_exp(1:53)), ...
    log(transfer_rev(1:53)), log(tax(1:53)), log(gdpQ(1:53)), ...
    spreads(1:53)*100, log(govc(1:53)), dum_eap(1:53), ...
    'VariableNames',{'T_t','Te_t','Tr_t','Tax_t','Y_t','RS_t','G_t','Dum_eap_t'}) ; 



%%
jcitest([tb_rs_predef.Te_t tb_rs_predef.Y_t  ]) 
% test if cointegration, r0=1 r1=0 means one cointegration relationship

md_trans_exp_predef_eap = fitlm(tb_rs_predef, 'Te_t ~ Y_t*Dum_eap_t ') ; 
disp(md_trans_exp_predef_eap )

res = md_trans_exp_predef_eap.Residuals.Raw ;

[h, pValue, stat, cValue, reg] = adftest(res) ; 
% test whether residuals are stationary, h=1 confirms stationary

%% test ECT error correction term

md_trans_exp_predef_eap = fitlm(tb_rs_predef, 'Te_t ~ Y_t*Dum_eap_t ') ; 
residuals = md_trans_exp_predef_eap.Residuals.Raw ; 
ECT = lagmatrix(residuals, 1); 
Delta_Te_t = diff(tb_rs_predef.Te_t);
Delta_Y_t = diff(tb_rs_predef.Y_t);

tb_ecm = table(Delta_Te_t, Delta_Y_t, ECT(2:end), tb_rs_predef.Dum_eap_t(2:end), ...
               'VariableNames', {'Delta_Te_t', 'Delta_Y_t', 'ECT', 'Dum_eap_t'});

% Fit the ECM
ecm_model = fitlm(tb_ecm, 'Delta_Te_t ~ Delta_Y_t*Dum_eap_t + ECT  ');
disp(ecm_model);


%%
jcitest([tb_rs_predef.Tr_t tb_rs_predef.Y_t  ]) 
md_trans_rev_predef_eap = fitlm(tb_rs_predef, 'Tr_t ~ Y_t*Dum_eap_t ') ; 
disp(md_trans_rev_predef_eap )

res = md_trans_rev_predef_eap.Residuals.Raw ;
[h, pValue, stat, cValue, reg] = adftest(res) ; 


%%
jcitest([tb_rs_predef.G_t tb_rs_predef.Y_t  ]) 

md_govc_predef_eap = fitlm(tb_rs_predef, 'G_t ~ Y_t + Dum_eap_t*Y_t ') ; 
disp(md_govc_predef_eap )

res = md_govc_predef_eap.Residuals.Raw;  % test residual stationarity
[h, pValue, stat, cValue, reg] = adftest(res) ; % h=1, stationary. 


%% How EAP affect taxes?
jcitest([tb_rs_predef.Tax_t tb_rs_predef.Y_t  ]) 

md_tax_predef_eap = fitlm(tb_rs_predef, 'Tax_t ~ Y_t + Dum_eap_t*Y_t ') ; 
disp(md_tax_predef_eap )

res = md_tax_predef_eap.Residuals.Raw;  % test residual stationarity
[h, pValue, stat, cValue, reg] = adftest(res) ; % h=1, stationary. 


%% How T_t affect spreads?

jcitest([tb_rs_predef.RS_t tb_rs_predef.Y_t  tb_rs_predef.Te_t ])

md_rs_predef = fitlm(tb_rs_predef, 'RS_t ~ Y_t  + Te_t ') ; 
disp(md_rs_predef )

res = md_rs_predef.Residuals.Raw;  % test residual stationarity
[h, pValue, stat, cValue, reg] = adftest(res) ;



md_rs_predef_eap = fitlm(tb_rs_predef, 'RS_t ~ Y_t + Te_t + Dum_eap_t:Te_t + Dum_eap_t:Y_t ') ; 
disp(md_rs_predef_eap )

res = md_rs_predef_eap.Residuals.Raw;  % test residual stationarity
[h, pValue, stat, cValue, reg] = adftest(res) ;





