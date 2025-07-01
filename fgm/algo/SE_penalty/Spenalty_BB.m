%% ****************************************************************
%   filename: Spenalty_BB
%% ****************************************************************

function [fobj,Xk,gapk,iter] = Spenalty_BB(Xk,retr,sqhfun,hfun,OPTIONS,K,n, r)

%%
if isfield(OPTIONS,'maxiter');      maxiter    = OPTIONS.maxiter;     end
if isfield(OPTIONS,'printyes');     printyes   = OPTIONS.printyes;    end
if isfield(OPTIONS,'gtol');         gtol       = OPTIONS.gtol;        end
if isfield(OPTIONS,'objtol');       objtol     = OPTIONS.objtol;      end
if isfield(OPTIONS,'objX');         objX       = OPTIONS.objX;        end

% [n,r] = size(K);

Kp = K; % Bp = B;

%% ***************************************************************

fk = objfun_QAP(Xk,K,n,r,objX);

gapk = sqhfun(Xk);

ratio = abs(fk)/gapk;

%% ****************** to estimate the initial rho *****************

if ratio >=1.0e+7

    rho = abs(fk)/(1.0e+8*gapk);
    
    gamma0 = 1e-8; 
    
elseif ratio >= 5.0e+5

    rho = abs(fk)/(1.0e+8*gapk);
    
    gamma0 = 1e-5;    
    
elseif ratio >= 5.0e+4

    rho = abs(fk)/(1.0e+6*gapk);
    
    gamma0 = 1e-3;
    
elseif ratio >= 5.0e+3

    rho = abs(fk)/(1.0e+5*gapk);
    
    gamma0 = 5e-3;
    
elseif ratio >= 5.0e+2

    rho = abs(fk)/(1.0e+5*gapk);
    
    gamma0 = 1e-2;
    
elseif ratio >= 5.0e+1

    rho = abs(fk)/(1.0e+5*gapk);
    
    gamma0 = 5e-2;
else
    rho = abs(fk)/(1.0e+4*gapk);
    
    gamma0 = 5e-2;
end

pars.max_gam = gamma0;

pars.min_gam = 1e-12;

%% ******************* parameters for dSNCG ************************

OPTIONS_SpBB.printyes = 0;

OPTIONS_SpBB.maxiter = 100;

OPTIONS_SpBB.Vtol = 5.0e-3;

OPTIONS_SpBB.Ftol = 1.0e-8;

OPTIONS_SpBB.objX = objX;
%% ***************************************************************

tstart = clock;

Fk_list = zeros(maxiter,1);

for iter = 1:maxiter
    
    pars.rho = rho;
    
    pars.gamma = gamma0;
    
%     [Xk,Fk,subiter] = PG_BB(Xk,retr,pars,OPTIONS_SpBB,K,n, r);
    [Xk,subiter,objfk] = PG_BB(Xk,retr,pars,OPTIONS_SpBB,K,n, r);
    gapk = hfun(Xk);
    
    Fk_list(iter) = objfk;
    
    ttime = etime(clock,tstart);
    
    if (printyes)
        
        fprintf('\n %3d    %3d       %3.2e      %5.4e     %3.2e     %2.1f',iter,subiter,gapk,Fk,rho,ttime);
        
    end
    
    %% *************** to check the stopping condition **************
    
    if (gapk<=gtol)
        
        fobj = objfun_QAP(Xk, Kp, n, r,objX);
        
        return;
        
    elseif (gapk<=5*gtol && iter>=10 && abs(objfk-Fk_list(iter-9))/(1+abs(objfk))<=objtol)
        
        fobj = objfun_QAP(Xk,Kp,n, r,objX);
        
        return;
    end
    
    rho = min(rho*1.1,1e+10);
    
    Xk = abs(Xk);
    
    Xk = myQR(Xk,r); 

end

fobj = objfun_QAP(Xk,Kp,n, r,objX);
