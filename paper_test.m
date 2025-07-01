%% *********** test the dPPA plus SNCG for QAP examples ************

clear;

close all;

restoredefaultpath;

addpath(genpath('QAPData'));

addpath(genpath('solvers'));

addpath(genpath('dSNCG'));

addpath(genpath('SE_penalty'));

addpath(genpath('OptM-master'));

% ***************** generate a QAP problem ***********************

load('qap_name');

load('best_val');

prob_No = 57

ss = strcat(qap_name(prob_No),'.mat');

load(ss);

bvalue = best_val(prob_No);

[n,r] = size(A);

Asnorm = svds(A,1);

Bsnorm = svds(B,1);
%% ******************* parameters for dSNCG *************************

OPTIONS1.maxiter = 2000;

OPTIONS1.printyes = 0;

OPTIONS1.gtol = 1e-3;

OPTIONS1.ftol = 1e-8;

%%
OPTIONS2.maxiter = 2000;

OPTIONS2.printyes = 1;

OPTIONS2.gtol = 1.0e-3;

OPTIONS2.objtol = 1.0e-8;

%%
opts.record = 0;
opts.tol = 1e-3;
opts.tolsub = 1e-3;
opts.omxitr = 100;
opts.mxitr = 100;
opts.xtol = 1e-5;
opts.gtol = 1e-5;
opts.ftol = 1e-8;
opts.tau = 1e-3;
opts.objX = 2;

%% ***************** to define the function ************************

sqhfun = @(x)sum(x(x<0).^2);      % sum(sum(max(0,-x).^2))

hfun = @(x)abs(sum(x(x<0)));      % sum(sum(max(0,-x)))

retr = @(x,v) retr_st(x,v,r,3);

printyes = 1;

ntest = 100;

gap_list1 = zeros(ntest,1);

orth_list1 = zeros(ntest,1);

error_list1 = zeros(ntest,1);

fobj_list1 = zeros(ntest,1);

time_list1 = zeros(ntest,1);

gap_list2 = zeros(ntest,1);

orth_list2 = zeros(ntest,1);

error_list2 = zeros(ntest,1);

fobj_list2 = zeros(ntest,1);

time_list2 = zeros(ntest,1);

gap_list3 = zeros(ntest,1);

orth_list3 = zeros(ntest,1);

error_list3 = zeros(ntest,1);

fobj_list3 = zeros(ntest,1);

time_list3 = zeros(ntest,1);

gap_list4 = zeros(ntest,1);

orth_list4 = zeros(ntest,1);

error_list4 = zeros(ntest,1);

fobj_list4 = zeros(ntest,1);

time_list4 = zeros(ntest,1);

for kk = 1:ntest
    
    randstate = 100+(kk-1)*9867;
    
    randn('state',double(randstate));
    
    rand('state',double(randstate));
    
    fprintf('\n  the kk = %2d problem',kk);
    
    state = rng;
    
    %% ********************** Initialization part ********************
    
    tempX0 = randn(n,r);
    
    X0 = orth(tempX0);
    
    rng(state);
    
%     tstart = clock;
    
%     [fobj1,gap,Xsol] = dSNCG_QAPBB(X0,retr,hfun,OPTIONS1,A,B,Asnorm,Bsnorm);
    
%     time_list1(kk) = etime(clock,tstart);
%     
%     orth_list1(kk) = norm(Xsol'*Xsol-eye(r),'fro')
%     
%     gap_list1(kk) = gap;
%     
%     fobj_list1(kk) = fobj1;
%     
%     error_list1(kk) = (fobj1-bvalue)/bvalue;
     
    %% *
    tstart = clock;
    
    state = rng;
    
    [fobj2,Xsol2] = Mpenalty_BB(X0,retr,sqhfun,hfun,OPTIONS2,A,B);
    
    orth_list2(kk) = norm(Xsol2'*Xsol2-eye(r),'fro');
    
    gap_list2(kk) = hfun(Xsol2);
    
    fobj_list2(kk) = fobj2;
    
    error_list2(kk) = (fobj2-bvalue)/bvalue;
    
    time_list2(kk) = etime(clock,tstart);
    
    rng(state);
    %%
%     tstart = clock;
%     
%     state = rng;
%     
%     [fobj3,Xsol3] = Spenalty_BB(X0,retr,sqhfun,hfun,OPTIONS2,A,B);
%     
%     orth_list3(kk) = norm(Xsol3'*Xsol3-eye(r),'fro');
%     
%     gap_list3(kk) = hfun(Xsol3);
%     
%     fobj_list3(kk) = fobj3;
%     
%     error_list3(kk) = (fobj3-bvalue)/bvalue;
%     
%     time_list3(kk) = etime(clock,tstart);
%     
%     rng(state);
    
%     tstart = clock;
%     
%     [Xsol4, out] = QAPORTHALEQQ_QUAD(X0,A,B,10,opts);
%     
%     orth_list4(kk) = norm(Xsol4'*Xsol4-eye(r),'fro');
%     
%     gap_list4(kk) = hfun(Xsol4);
%     
%     fobj_list4(kk) = objfun_QAP(Xsol4,A,B);
%     
%     error_list4(kk) = (fobj_list4(kk)-best_val(prob_No))/best_val(prob_No);
%     
%     time_list4(kk) = etime(clock,tstart);

end
max(orth_list1)
SNCG_gap = max(gap_list1)
SNCG_maxerr = max(error_list1)
SNCG_minerr = min(error_list1)
SNCG_mederr = mean(error_list1)
SNCG_time = mean(time_list1)
max(orth_list2)
NEP_gap = max(gap_list2)
NEP_maxerr = max(error_list2)
NEP_minerr = min(error_list2)
NEP_mederr = mean(error_list2)
NEP_time = mean(time_list2)
max(orth_list3)
sqEP_gap = max(gap_list3)
sqEP_maxerr = max(error_list3)
sqEP_minerr = min(error_list3)
sqEP_mederr = mean(error_list3)
sqEP_time = mean(time_list3)
max(orth_list4)
ALM_gap = max(gap_list4)
ALM_maxerr = max(error_list4)
ALM_minerr = min(error_list4)
ALM_mederr = mean(error_list4)
ALM_time = mean(time_list4)
