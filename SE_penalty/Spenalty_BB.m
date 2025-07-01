%% ****************************************************************
%   filename: Spenalty_BB
%% ****************************************************************

function [fobj,Xk,gapk,iter] = Spenalty_BB(Xk,retr,sqhfun,hfun,OPTIONS,A,B)

%%
if isfield(OPTIONS,'maxiter');      maxiter    = OPTIONS.maxiter;     end
if isfield(OPTIONS,'printyes');     printyes   = OPTIONS.printyes;    end
if isfield(OPTIONS,'gtol');         gtol       = OPTIONS.gtol;        end
if isfield(OPTIONS,'objtol');       objtol     = OPTIONS.objtol;      end

[n,r] = size(A);

Ap = A;  Bp = B;

%% ***************************************************************

fk = objfun_QAP(Xk,A,B);

gapk = sqhfun(Xk);

abs(fk)/gapk

%% ****************** to estimate the initial rho *****************

if abs(fk)/gapk>=1.0e+7

    rho = abs(fk)/(1.0e+8*gapk)
    
    gamma0 = 1e-8; 
    
elseif abs(fk)/gapk>=5.0e+5

    rho = abs(fk)/(1.0e+8*gapk)
    
    gamma0 = 1e-5;    
    
elseif abs(fk)/gapk>=5.0e+4

    rho = abs(fk)/(1.0e+6*gapk)
    
    gamma0 = 1e-3;
    
elseif abs(fk)/gapk>=5.0e+3

    rho = abs(fk)/(1.0e+5*gapk);
    
    gamma0 = 5e-3;
    
elseif abs(fk)/gapk>=5.0e+2

    rho = abs(fk)/(1.0e+5*gapk)
    
    gamma0 = 1e-2;
    
elseif abs(fk)/gapk>=5.0e+1

    rho = abs(fk)/(1.0e+5*gapk)
    
    gamma0 = 5e-2;
else
    rho = abs(fk)/(1.0e+4*gapk)
    
    gamma0 = 5e-2;
end

pars.max_gam = gamma0;

pars.min_gam = 1e-12;

%% ******************* parameters for dSNCG ************************

OPTIONS_SpBB.printyes = 0;

OPTIONS_SpBB.maxiter = 100;

OPTIONS_SpBB.Vtol = 5.0e-3;

OPTIONS_SpBB.Ftol = 1.0e-8;

%% ***************************************************************

tstart = clock;

fk_list = zeros(maxiter,1);

for iter = 1:maxiter
    
    pars.rho = rho;
    
    pars.gamma = gamma0;
    
    [Xk,subiter,objfk] = PG_BB(Xk,retr,pars,OPTIONS_SpBB,A,B);
    
    gapk = hfun(Xk);
    
    fk_list(iter) = objfk; 
    
    ttime = etime(clock,tstart);
    
    if (printyes)
        
        fprintf('\n %3d    %3d       %3.2e      %5.4e     %3.2e     %2.1f',iter,subiter,gapk,objfk,rho,ttime);
        
    end
    
  %% ***************** to check the stopping condition **************
    
    if (gapk<=gtol)
        
        fobj = objfun_QAP(Xk,Ap,Bp);
        
        return;
        
    elseif (gapk<=5*gtol && iter>=10 && abs(objfk-fk_list(iter-9))/(1+abs(objfk))<=objtol)
        
        fobj = objfun_QAP(Xk,Ap,Bp);
        
        return;
    end
    
    rho = min(rho*1.1,1e+10);
    
    Xk = abs(Xk)+1.0e-2*randn(n,r);
    
    Xk = myQR(Xk,r); 

end

fobj = objfun_QAP(Xk,Ap,Bp);
