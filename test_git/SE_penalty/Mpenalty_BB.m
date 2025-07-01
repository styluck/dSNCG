%% ****************************************************************
%   filename: Mpenalty_BB
%% ****************************************************************

function [fobj,Xk,gapk,iter] = Mpenalty_BB(Xk,OPTIONS,varargin)
%%

if isfield(OPTIONS,'maxiter');    maxiter    = OPTIONS.maxiter;     end
if isfield(OPTIONS,'printyes');   printyes   = OPTIONS.printyes;    end
if isfield(OPTIONS,'gtol');       gtol       = OPTIONS.gtol;        end
if isfield(OPTIONS,'ftol');       objtol     = OPTIONS.objtol;      end

[n,r] = size(Xk);

%% ***************************************************************
retr = @(x,v) retr_st(x,v,r,3);
sqhfun = @(x)sum(x(x<0).^2);      % sum(sum(max(0,-x).^2))
hfun = @(x)abs(sum(x(x<0)));      % sum(sum(max(0,-x)))
fk = objfun(Xk,varargin{1});

gapk = sqhfun(Xk); 

%abs(fk)/gapk

%% ****************** to estimate the initial rho *****************

if abs(fk)/gapk>=1.0e+7

    rho = abs(fk)/(1.0e+8*gapk);
    
    gamma0 = 1e-6; 
    
elseif abs(fk)/gapk>=5.0e+5

    rho = abs(fk)/(1.0e+8*gapk);
    
    gamma0 = 1e-5;    
    
elseif abs(fk)/gapk>=5.0e+4
     
    rho = abs(fk)/(1.0e+6*gapk);
    
    gamma0 = 1e-3;
    
elseif abs(fk)/gapk>=5.0e+3
    
    rho = abs(fk)/(1.0e+5*gapk);
    
    gamma0 = 5e-3;
    
elseif abs(fk)/gapk>=5.0e+2

    rho = abs(fk)/(1.0e+5*gapk);
    
    gamma0 = 1e-2;
    
elseif abs(fk)/gapk>=5.0e+1

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

OPTIONS_SpBB.maxiter = 1000;

OPTIONS_SpBB.Vtol = 5.0e-3;

OPTIONS_SpBB.Ftol = 1.0e-8;

%% ***************************************************************

tstart = clock;

Fk_list = zeros(maxiter,1);

for iter = 1:maxiter
    
    pars.rho = rho;
    
    pars.gamma = gamma0;
    
    [Xk,Fk,subiter] = MPG_BB(Xk,retr,pars,OPTIONS_SpBB,varargin{:});
    
    gapk = hfun(Xk);
    
    Fk_list(iter) = Fk;
    
    ttime = etime(clock,tstart);
    
    if (printyes)
        
        fprintf('\n %3d    %3d       %3.2e      %5.4e     %3.2e     %2.1f',iter,subiter,gapk,Fk,rho,ttime);
        
    end
    
    %% *************** to check the stopping condition **************
    
    if (gapk<=gtol)

        fobj = objfun(Xk,varargin{1});
        
        return;
        
    elseif (gapk<=5*gtol && iter>=10 && abs(Fk-Fk_list(iter-9))/(1+abs(Fk))<=objtol)
    
        fobj = objfun(Xk,varargin{1});
        
        return;
    end
    
    rho = min(rho*1.1,1e+10);
    
    Xk = abs(Xk)+1.0e-2*randn(n,r);
    
    Xk = myQR(Xk,r); 

end

fobj = objfun(Xk,varargin{1});
