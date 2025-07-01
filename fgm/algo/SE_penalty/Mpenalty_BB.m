%% ****************************************************************
%   filename: Mpenalty_BB
%% ****************************************************************

function [fobj,Xk,gapk,iter] = Mpenalty_BB(Xk,retr,hfun,OPTIONS,K,n,r)

%%
if isfield(OPTIONS,'maxiter');      maxiter    = OPTIONS.maxiter;     end
if isfield(OPTIONS,'printyes');     printyes   = OPTIONS.printyes;    end
if isfield(OPTIONS,'gtol');         gtol       = OPTIONS.gtol;        end
if isfield(OPTIONS,'objtol');       objtol     = OPTIONS.objtol;      end
if isfield(OPTIONS,'objX');         objX       = OPTIONS.objX;        end

%% ***************************************************************

fk = objfun_QAP(Xk,K,n,r,objX);

gapk = hfun(Xk);

ratio = abs(fk)/gapk;

%% ****************** to estimate the initial rho *****************

if ratio>=1.0e+7

    rho = abs(fk)/(1.0e+8*gapk);
    
    gamma0 = 1e-6; 
    
elseif ratio>=5.0e+5

    rho = abs(fk)/(1.0e+7*gapk);
    
    gamma0 = 1e-5;    
    
elseif ratio>=5.0e+4
     
    rho = abs(fk)/(1.0e+5*gapk);
    
    gamma0 = 1e-3;
    
elseif ratio>=5.0e+3
    
    rho = abs(fk)/(1.0e+5*gapk);
    
    gamma0 = 5e-3;
    
elseif ratio>=5.0e+2

    rho = abs(fk)/(1.0e+5*gapk);
    
    gamma0 = 1e-2;
    
else
    rho = abs(fk)/(1.0e+4*gapk);
    
    gamma0 = 5e-2;
%     
%     rho = abs(fk)/(1.0e+6*gapk)
%     
%     gamma0 = 1e+4;
end

pars.max_gam = gamma0;

pars.min_gam = 1e-12;

%% ******************* parameters for dSNCG ************************

OPTIONS_MpBB.printyes = 0;

OPTIONS_MpBB.maxiter = 1000;

OPTIONS_MpBB.Vtol = 5.0e-3;

OPTIONS_MpBB.Ftol = 1.0e-8;

OPTIONS_MpBB.objX = objX;
%% ***************************************************************

tstart = clock;

fk_list = zeros(maxiter,1);

for iter = 1:maxiter
    
    pars.rho = rho;
    
    pars.gamma = gamma0;
    
    [Xk,subiter,objfk] = MPG_BB(Xk,retr,pars,OPTIONS_MpBB,K,n,r);
    
    gapk = hfun(Xk);
    
    fk_list(iter) = objfk;
    
    ttime = etime(clock,tstart);
    
    if (printyes)
        
        fprintf('\n %3d    %3d       %3.2e      %5.4e     %3.2e     %2.1f',iter,subiter,gapk,objfk,rho,ttime);
        
    end
    
    %% *************** to check the stopping condition **************
    
    if (gapk<=gtol)

        fobj = objfun_QAP(Xk,K,n,r,objX);
        
        return;
        
    elseif (gapk<=5*gtol && iter>=10 && abs(objfk-fk_list(iter-9))/(1+abs(objfk))<=objtol)
    
        fobj = objfun_QAP(Xk,K,n,r,objX);
        
        return;
    end
    
    rho = min(rho*1.1,1e+10);
    
    Xk = abs(Xk);%+min(1e-2,gapk)*randn(n,r);
    
    Xk = myQR(Xk,r); 
    
    OPTIONS_MpBB.Vtol = max(0.99*OPTIONS_MpBB.Vtol,1.0e-3);

end

fobj = objfun_QAP(Xk,K,n,r,objX);
