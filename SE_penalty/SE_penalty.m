%% ****************************************************************
%   filename: SE_penalty
%% ****************************************************************
% the minimization of smooth exact penalty over the Stiefel manifold

function [fobj,Xk,iter] = SE_penalty(Xk,sqhfun,hfun,OPTIONS,A,B)

%%
if isfield(OPTIONS,'maxiter');          maxiter    = OPTIONS.maxiter;   end
if isfield(OPTIONS,'printyes');         printyes   = OPTIONS.printyes;  end
if isfield(OPTIONS,'gtol');             gtol       = OPTIONS.gtol;      end
if isfield(OPTIONS,'ftol');             ftol       = OPTIONS.ftol;      end

[n,r] = size(A);

opts.record = 0;
opts.mxitr = 100;
opts.xtol = 1e-5;
opts.gtol = 1e-5;
opts.ftol = 1e-8;
opts.eta  = 0.1;
opts.tau = 1e-3;

fk = objfun_QAP(Xk,A,B);

gapk = sqhfun(Xk);

if abs(fk)/gapk>=1.0e+4
    
    rho = abs(fk)/(1.0e+6*gapk);   %% this is the best
    
elseif abs(fk)/gapk>=1.0e+3
    
    rho = abs(fk)/(1.0e+3*gapk);

elseif abs(fk)/gapk>=1.0e+2
    
    rho = abs(fk)/(5.0e+2*gapk);    
    
elseif abs(fk)/gapk>=1.0e+1
    
    rho = abs(fk)/(1.0e+2*gapk);
else
    rho = abs(fk)/(10*gapk);
end

tstart = clock;

Fk_list = zeros(maxiter,1);

for iter = 1:maxiter
    
    gapk_old = gapk;
    
    [Xk,out]= OptStiefelGBB(Xk, @Epen_fgrad,opts,A,B,rho);
    
    Fk = out.fval;
    
    Fk_list(iter) = Fk;
    
    gapk = hfun(Xk);   % do not change it into hfun1
    
    ratio = gapk/gapk_old;
         
    iter_FOptM = out.itr;    
    
    ttime = etime(clock,tstart);
    
    if (printyes)
        
        fk = objfun_QAP(Xk,A,B);
        
        fprintf('\n %3d    %3d       %3.2e      %6.5e     %3.2e     %2.1f   %3.2f',iter,iter_FOptM,gapk,fk,rho,ttime,ratio);
        
    end
   
    if gapk<gtol
        
       break;
        
    elseif (gapk<=5*gtol && iter>=10 && abs(Fk-Fk_list(iter-9))/(1+abs(Fk))<=ftol)
       
       break;        
    end
 %% ********************* update the parameter rho ****************    
    if ratio==1        
        rho = min(rho*1.2,1e+10);
    else
        rho = min(rho*1.1,1e+10);
    end
    
    Xk = abs(Xk)+1.0e-2*randn(n,r);    %% This is important !!
    
end
if (printyes)
    fobj = fk;
else
    fobj = objfun_QAP(Xk,A,B);
end
