%% ****************************************************************
%   filename: dSNCG_QAPBB
%% ****************************************************************

function [fobj,gapk,Xk,iter] = dSNCG_QAPBB(Xk,OPTIONS,varargin)

%%
if isfield(OPTIONS,'maxiter');      maxiter    = OPTIONS.maxiter;     end
if isfield(OPTIONS,'printyes');     printyes   = OPTIONS.printyes;    end
if isfield(OPTIONS,'gtol');         gtol       = OPTIONS.gtol;        end
if isfield(OPTIONS,'ftol');         ftol       = OPTIONS.ftol;        end

[n,r]=size(Xk);

max_norm = svds(varargin{1},1);

%% ***************************************************************
retr = @(x,v) retr_st(x,v,r,3);
hfun = @(x)abs(sum(x(x<0)));      % sum(sum(max(0,-x)))

[fk,gradk] = objfun(Xk,varargin{1});

gapk = hfun(Xk);

if abs(fk)/gapk<5.0e+1
    
    if max_norm<=10 || abs(fk)/gapk<=10
        varargin{1} = varargin{1}*1.0e+3;
    else
        varargin{1} = varargin{1}*1.0e+1;
    end
    
    [fk,gradk] = objfun(Xk,varargin{1});
    
     gapk = hfun(Xk);
end

% abs(fk)/gapk

%% ****************** to estimate the initial rho *****************

if abs(fk)/gapk>=1.0e+7

    rho = abs(fk)/(1.0e+8*gapk);
    
    gamma0 = 1e-7;
    
elseif abs(fk)/gapk>=1.0e+6

    rho = abs(fk)/(1.0e+7*gapk);
    
    gamma0 = 1e-6;   % do not increase it !!
    
elseif abs(fk)/gapk>=1.0e+5

    rho = abs(fk)/(1.0e+6*gapk);
    
    gamma0 = 1e-5;
    
elseif abs(fk)/gapk>=1.0e+4

    rho = abs(fk)/(1.0e+5*gapk);
    
    gamma0 = 5e-5;
    
elseif abs(fk)/gapk>=5.0e+2

    rho = abs(fk)/(1.0e+4*gapk);
    
    gamma0 = 1e-3;
    
else 
% elseif abs(fk)/gapk>=5.0e+1
    rho = abs(fk)/(1.0e+2*gapk);
    
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
    
    [U,Xk,fk,gapk,subiter] = dSNCG_BB(U,Xk,gradk,fk,hfun,retr,pars,OPTIONS_dSNCG,varargin{1},r);

    fk_list(iter) = fk;
    
    ttime = etime(clock,tstart);
    
    if (printyes)
        
        fprintf('\n %3d    %3d       %3.2e      %5.4e     %3.2e     %2.1f',iter,subiter,gapk,fk,rho,ttime);
        
    end
    
    %% *************** to check the stopping condition **************
    
    if (gapk<=gtol)

        fobj = objfun(Xk,varargin{1});
        
        return;
        
    elseif (gapk<=5*gtol && iter>=10 && abs(fk-fk_list(iter-9))/(1+abs(fk))<=ftol)
        
        fobj = objfun(Xk,varargin{1});
        
        return;
    end

    rho = min(rho*1.1,1e+10);

    Xk = myQR(abs(Xk),r);

    [fk,gradk] = objfun(Xk,varargin{1});
end

fobj = objfun(Xk,varargin{1});