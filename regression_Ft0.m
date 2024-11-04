clear;  clc
format compact

simu_longpath

highstress = 0.95; % 0.965
lowstress = 0.05;  % 0.035

BadD = Bad+D;
FE = FE*mean(F)*100 ;  % in percentage points
Z = Z*mean(Z)*100 ;
B = B*mean(B)*100 ;
PM = PM*100; % in basis points

%% 1. select data in good financial status
n = length(PM);
pm = nan(n,1);
fe = pm; zz = pm;  bb = pm; ind = pm ;

for i = 3:n
    if  BadD(i-1)~=1 && BadD(i-2)~=1 && BadD(i)~=1
       pm(i) = PM(i);
       fe(i) = FE(i);
       zz(i) = Z(i);
       bb(i) = B(i);
    end
end

pm = pm(~isnan(pm));  % in basis points
fe = fe(~isnan(fe));  % in percentage
zz = zz(~isnan(zz));
bb = bb(~isnan(bb));

[cdf,cdf_pm] = ecdf(pm);

figure
plot(cdf_pm/100,cdf,linewidth=2)
xlabel('Spreads (%)')
ylabel('Cumulative Density')
ylim([-0.01 1.01])
xlim([0 33])

tb1 = table(pm,fe,zz,bb,'VariableNames',{'Spreads','Fe_t','Z_t','B_t'});
md1 = fitlm(tb1,'Spreads ~ Fe_t + Z_t + B_t') ;  
disp('Unconditional')
disp(md1)


%% unconditional, cuts only
n = length(PM);
pm = nan(n,1);
fe = pm; zz = pm; bb = pm;

for i = 3:n
    if  BadD(i-1)~=1 && BadD(i-2)~=1 && BadD(i)~=1 && FE(i) < 0
       pm(i) = PM(i);
       fe(i) = FE(i);
       zz(i) = Z(i);
       bb(i) = B(i);
    end
end

pm = pm(~isnan(pm));  % in basis points
fe = fe(~isnan(fe));  % in percentage
zz = zz(~isnan(zz));
bb = bb(~isnan(bb));

tb11 = table(pm,fe,zz,bb,'VariableNames',{'Spreads','Fen_t','Z_t','B_t'});
md11 = fitlm(tb11,'Spreads ~ Fen_t + Z_t + B_t') ;  
disp('Uncondtional, cuts only')
disp(md11)


%% unconditional, rise only
n = length(PM);
pm = nan(n,1);
fe = pm; zz = pm; bb = pm; 

for i = 3:n
    if  BadD(i-1)~=1 && BadD(i-2)~=1 && BadD(i)~=1 && FE(i) > 0
       pm(i) = PM(i);
       fe(i) = FE(i);
       zz(i) = Z(i);
       bb(i) = B(i);
    end
end

pm = pm(~isnan(pm));  % in basis points
fe = fe(~isnan(fe));  % in percentage
zz = zz(~isnan(zz));
bb = bb(~isnan(bb));

tb12 = table(pm,fe,zz,bb,'VariableNames',{'Spreads','Fep_t','Z_t','B_t'});
md12 = fitlm(tb12,'Spreads ~ Fep_t + Z_t + B_t') ;  
disp('Unconditional, rises only')
disp(md12)


%%  2. estimate the case of high financial stress

lower_bound = min(cdf_pm(cdf >= highstress));  % define the level of "high stress"

% write a loop to pick data whose cdf is bigger/smaller than a threshold
n = length(PM);
pm1 = nan(n,1);
fe1 = pm1; zz1 = pm1; bb1 = pm1;

for i = 3:n
    if  PM(i-1) >= lower_bound && BadD(i-1)~=1 && BadD(i-2)~=1 && BadD(i)~=1
       pm1(i) = PM(i);
       fe1(i) = FE(i);
       zz1(i) = Z(i);
       bb1(i) = B(i);
    end
end

pm1 = pm1(~isnan(pm1));  % in basis points
fe1 = fe1(~isnan(fe1));  % in percentage
zz1 = zz1(~isnan(zz1));
bb1 = bb1(~isnan(bb1));

tb2 = table(pm1,fe1,zz1,bb1,'VariableNames',{'Spreads','Fe_t','Z_t','B_t'});
md2 = fitlm(tb2,'Spreads ~ Fe_t + Z_t + B_t') ;  
disp('High stress')
disp(md2)


%% 3. estimate the case of low financial stress

upper_bound = max(cdf_pm(cdf <= lowstress));  % define the level of "low stress"

% write a loop to pick data whose cdf is bigger/smaller than a threshold
n = length(pm);
pm1 = nan(n,1);
fe1 = pm1; zz1 = pm1; bb1 = pm1;

for i = 3:n
    if  PM(i-1) <= upper_bound && BadD(i-1)~=1 && BadD(i-2)~=1 && BadD(i)~=1
       pm1(i) = PM(i);
       fe1(i) = FE(i);
       zz1(i) = Z(i);
       bb1(i) = B(i);
    end
end

pm1 = pm1(~isnan(pm1));  % in basis points
fe1 = fe1(~isnan(fe1));  % in percentage
zz1 = zz1(~isnan(zz1));
bb1 = bb1(~isnan(bb1));

tb3 = table(pm1,fe1,zz1,bb1,'VariableNames',{'Spreads','Fe_t','Z_t','B_t'});
md3 = fitlm(tb3,'Spreads ~ Fe_t + Z_t + B_t') ;  
disp('Low stress, rises only')
disp(md3)


%% 4. Split G_t with increases and cuts, with high financial stress

n = length(pm);
pm_p = nan(n,1);    pm_n = nan(n,1);
fe_p = nan(n,1);    fe_n = nan(n,1);
zz_p = nan(n,1);    zz_n = nan(n,1);
bb_p = nan(n,1);    bb_n = nan(n,1);    

for i = 3:n
    if  PM(i-1) >= lower_bound && BadD(i-1)~=1 && BadD(i-2)~=1  && FE(i) < 0 && BadD(i)~=1 
        zz_n(i) = Z(i);
        bb_n(i) = B(i);
        pm_n(i) = PM(i);
        fe_n(i) = FE(i);
    end
    if PM(i-1) >= lower_bound && BadD(i-1)~=1 && BadD(i-2)~=1  && FE(i) > 0 && BadD(i)~=1 
        zz_p(i) = Z(i);
        bb_p(i) = B(i);
        pm_p(i) = PM(i);
        fe_p(i) = FE(i);
    end
end

pm_p = pm_p(~isnan(pm_p));  % in basis points
pm_n = pm_n(~isnan(pm_n));  % in basis points
fe_p = fe_p(~isnan(fe_p));   fe_n = fe_n(~isnan(fe_n)); % in percentage
zz_p = zz_p(~isnan(zz_p));   zz_n = zz_n(~isnan(zz_n));
bb_p = bb_p(~isnan(bb_p));   bb_n = bb_n(~isnan(bb_n));

tb4_n = table(pm_n,fe_n,zz_n,bb_n,'VariableNames',{'RSn_t','Fen_t','Zn_t','Bn_t'});
md4_n = fitlm(tb4_n,'RSn_t ~ Fen_t + Zn_t + Bn_t') ;  
disp('High stress, cuts only')
disp(md4_n)



tb4_p = table(pm_p,fe_p,zz_p,bb_p,'VariableNames',{'RSp_t','Fep_t','Zp_t','Bp_t'});
md4_p = fitlm(tb4_p,'RSp_t ~ Fep_t + Zp_t + Bp_t') ;  
disp('High stress, rises only')
disp(md4_p)



%% 5. Cuts/rise only, with low financial stress

n = length(pm);
pm_p = nan(n,1);    pm_n = nan(n,1);
fe_p = nan(n,1);    fe_n = nan(n,1);
zz_p = nan(n,1);    zz_n = nan(n,1);
bb_p = nan(n,1);    bb_n = nan(n,1);    

for i = 3:n
    if  PM(i-1) <= upper_bound && BadD(i-1)~=1 && BadD(i-2)~=1 && FE(i) < 0 && BadD(i)~=1 
        zz_n(i) = Z(i);
        bb_n(i) = B(i);
        pm_n(i) = PM(i);
        fe_n(i) = FE(i);
    end
    if pm(i-1) <= upper_bound && BadD(i-1)~=1 && BadD(i-2)~=1 && FE(i) > 0 && BadD(i)~=1 
        zz_p(i) = Z(i);
        bb_p(i) = B(i);
        pm_p(i) = PM(i);
        fe_p(i) = FE(i);
    end
end

pm_p = pm_p(~isnan(pm_p));  % in basis points
pm_n = pm_n(~isnan(pm_n));  % in basis points
fe_p = fe_p(~isnan(fe_p));   fe_n = fe_n(~isnan(fe_n)); % in percentage
zz_p = zz_p(~isnan(zz_p));   zz_n = zz_n(~isnan(zz_n));
bb_p = bb_p(~isnan(bb_p));   bb_n = bb_n(~isnan(bb_n));

tb5_n = table(pm_n,fe_n,zz_n,bb_n,'VariableNames',{'RSn_t','Fen_t','Zn_t','Bn_t'});
md5_n = fitlm(tb5_n,'RSn_t ~ Fen_t + Zn_t + Bn_t') ;  
disp('Low stress, cuts only')
disp(md5_n)

[~, hacSE_md5n, hacCoeff_md5n] = hac(md5_n,Type="HC",Weights="HC0", Display="off") ;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tb5_p = table(pm_p,fe_p,zz_p,bb_p,'VariableNames',{'RSp_t','Fep_t','Zp_t','Bp_t'});
md5_p = fitlm(tb5_p,'RSp_t ~ Fep_t + Zp_t + Bp_t') ;  
disp('Low stress, rises only')
disp(md5_p)


md1_coeff = md1.Coefficients;  md11_coeff = md11.Coefficients;  md12_coeff = md12.Coefficients;
md2_coeff = md2.Coefficients;  md3_coeff = md3.Coefficients;    md4n_coeff = md4_n.Coefficients;
md4p_coeff = md4_p.Coefficients;  md5n_coeff = md5_n.Coefficients;
md5p_coeff = md5_p.Coefficients;




