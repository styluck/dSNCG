%% *****************************************************************
%  filename: lsearch
%% *****************************************************************
% function alp = line_search(Xk, Vk, f, fold, retr, params, lsopts)
% line search method with the following condition:
% while:
%           f(retr_{Xk}(aVk)) < f(X_k) - (a/2t)||Vk||_F^2
% inputs:
%   X_k: iteration of the current loop
%   V_k: current direction
%   retr: the analytic function of retraction
%   lsopts: struct with:
%           printyes: print loop, yes=1(default) no=0
%           t: param for Armijo condition, default 1e-3 (0~0.5)
%           tau: param for Armijo condition. default .9
%           tol: tolerance of the max iteration (default = 1e3)

function [Xnew,fnew,Cnew,gap_new,lsearch_flag,iter] = lsearch(Xk,Vk,sqnormVk,Fold,hfun,retr,lsopts,rho,A,B)

if isfield(lsopts,'maxiter');          maxiter    = lsopts.maxiter;    end
if isfield(lsopts,'printyes');         printyes   = lsopts.printyes;   end
if isfield(lsopts,'tstep');            tstep      = lsopts.tstep;      end
if isfield(lsopts,'beta');             beta       = lsopts.beta;       end


sqnormVk = sqnormVk/(2*tstep);

alpha = 1;  lsearch_flag = 0;

Xnew = retr(Xk,alpha*Vk);

Xq = Xnew.^2;  

AXB = A'*Xq*B;

fnew = sum(dot(Xq,AXB));

gap_new = hfun(Xnew);

Fnew = fnew+rho*gap_new;

for iter=1:maxiter
    
    if Fnew<= Fold - alpha*sqnormVk
        
        Cnew = 2*(AXB+(A*Xq*B')).*Xnew;
       
        Cnew = ProjXZ(Xnew,Cnew);
                       
        return;        
    end
    
    alpha = beta*alpha;
    
    lsearch_flag = 1;
    
    Xnew = retr(Xk,alpha*Vk);
    
    Xq = Xnew.^2;   
    
    AXB = A'*Xq*B;
    
    fnew = sum(dot(Xq,AXB));
    
    gap_new = hfun(Xnew);
    
    Fnew = fnew+rho*gap_new;
    
    if (printyes)
        
        fprintf('\n %3d     %3.2e      %3.2e       %2.1f',iter,alpha,Fnew,ttime);
        
    end
    
end

Cnew = 2*(AXB+(A*Xq*B')).*Xnew;

Cnew = ProjXZ(Xnew,Cnew);   % it is necessary!!
    

