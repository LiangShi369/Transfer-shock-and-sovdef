clear;  clc

format compact

% Greece data from 1995 to the beginning of 2019

load data.mat
%% Calculate statistics

%% A. Estimate fiscal fule: G_t = psi0 + psi1*Y_t + psi2*D_t

x = isnan(debtQ);

g = gQsa(~x);
y = gdpQ(~x);
d = debtQ(~x) ; 

g = recover(g,x);
y = recover(y,x);
d = recover(d,x);

%% model 1, with debt ajustment
jcitest([g y d])  

tb1 = table(g,y,d,'VariableNames',{'G_t','Y_t','D_t'}) ; 
md1 = fitlm(tb1, 'G_t~Y_t+D_t' );  % md1 = fitlm(X,gyQ,'RobustOpts','on')  
disp(md1)

res = md1.Residuals.Raw;
[h, pValue, stat, cValue, reg] = adftest(res ) ;


%% income tax rule
jcitest([tyQ/100 tfpQ/mean(tfpQ)]) 

tb5 = table(tyQ/100, tfpQ/mean(tfpQ),'VariableNames',{'TAU_t','TFP_t'}) ;
md5 = fitlm(tb5, 'TAU_t~TFP_t');  % md2 = fitlm(y_cyc,gyQ,'RobustOpts','on') 
disp(md5)

res = md5.Residuals.Raw;
[h, pValue, stat, cValue, reg] = adftest(res ) ;

