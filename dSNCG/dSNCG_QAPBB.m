%% ****************************************************************
%   filename: dSNCG_QAPBB
%% ****************************************************************

function [fobj,gapk,Xk,iter] = dSNCG_QAPBB(Xk,retr,hfun,OPTIONS,A,B,Asnorm,Bsnorm)

%%
if isfield(OPTIONS,'maxiter');      maxiter    = OPTIONS.maxiter;     end
if isfield(OPTIONS,'printyes');     printyes   = OPTIONS.printyes;    end
if isfield(OPTIONS,'gtol');         gtol       = OPTIONS.gtol;        end
if isfield(OPTIONS,'ftol');         ftol       = OPTIONS.ftol;        end

[n,r]=size(A);

Ap = A;  Bp = B; 

max_norm = max(Asnorm,Bsnorm);

%% ***************************************************************

[fk,gradk] = objfun_QAP(Xk,A,B);

gapk = hfun(Xk);

if abs(fk)/gapk<5.0e+1
    
    if max_norm<=10 || abs(fk)/gapk<=10
        if Asnorm>Bsnorm
            B = B*1.0e+3;
        else
            A = A*1.0e+3;
        end
    else
        if Asnorm>Bsnorm
            B = B*1.0e+1;
        else
            A = A*1.0e+1;
        end
    end
    
    [fk,gradk] = objfun_QAP(Xk,A,B);
    
     gapk = hfun(Xk);
end

abs(fk)/gapk

%% ****************** to estimate the initial rho *****************

if abs(fk)/gapk>=1.0e+7

    rho = abs(fk)/(1.0e+8*gapk)
    
    gamma0 = 1e-7;
    
elseif abs(fk)/gapk>=1.0e+6

    rho = abs(fk)/(1.0e+7*gapk)
    
    gamma0 = 1e-6;   % do not increase it !!
    
elseif abs(fk)/gapk>=1.0e+5

    rho = abs(fk)/(1.0e+6*gapk)
    
    gamma0 = 1e-5;
    
elseif abs(fk)/gapk>=1.0e+4

    rho = abs(fk)/(1.0e+5*gapk);
    
    gamma0 = 5e-5;
    
elseif abs(fk)/gapk>=5.0e+2

    rho = abs(fk)/(1.0e+4*gapk)
    
    gamma0 = 1e-3;
    
elseif abs(fk)/gapk>=5.0e+1

    rho = abs(fk)/(1.0e+2*gapk)
    
    gamma0 = 1e-2;
end

pars.max_gam = gamma0;

pars.min_gam = 1e-12;  %% do not increase it !!

%% ******************* parameters for dSNCG ************************

OPTIONS_dSNCG.printyes = 0;

OPTIONS_dSNCG.maxiter = 20;

OPTIONS_dSNCG.Vtol = 5.0e-3;

OPTIONS_dSNCG.Ftol = 1.0e-8;

%% ***************************************************************

tstart = clock;

fk_list = zeros(maxiter,1);

U = zeros(r);

for iter = 1:maxiter
    
    pars.rho = rho;
    
    pars.gamma = gamma0;
    
    [U,Xk,fk,gapk,subiter] = dSNCG_BB(U,Xk,gradk,fk,hfun,retr,pars,OPTIONS_dSNCG,A,B,r);

    fk_list(iter) = fk;
    
    ttime = etime(clock,tstart);
    
    if (printyes)
        
        fprintf('\n %3d    %3d       %3.2e      %5.4e     %3.2e     %2.1f',iter,subiter,gapk,fk,rho,ttime);
        
    end
    
    %% *************** to check the stopping condition **************
    
    if (gapk<=gtol)

        fobj = objfun_QAP(Xk,Ap,Bp);
        
        return;
        
    elseif (gapk<=5*gtol && iter>=10 && abs(fk-fk_list(iter-9))/(1+abs(fk))<=ftol)
        
        fobj = objfun_QAP(Xk,Ap,Bp);
        
        return;
    end

    rho = min(rho*1.1,1e+10);

    Xk = myQR(abs(Xk),r);

    [fk,gradk] = objfun_QAP(Xk,A,B);
end

fobj = objfun_QAP(Xk,Ap,Bp);