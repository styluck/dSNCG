%% *******************************************************************
%  filename: SNewton
%
%% *******************************************************************
%% Semismooth Newton Method for solving the system 
%
% nabla Psi(u)=0
%
% with Psi(u)=-AXk+AProx_{gamma*rho}h(Zk+gamma*AtU)=0 
% 
%% ***************************************************************
%%  

function [U,AtU,proxh,iter]= SNewton(U,AtU,pars,OPTIONS,Zk,Xk,AXk)

%% *************** Initialization of parameter *********************

if isfield(OPTIONS,'tol');           tol_SNCG   = OPTIONS.tol;         end
if isfield(OPTIONS,'printyes');      printyes   = OPTIONS.printyes;    end
if isfield(OPTIONS,'maxiter');       maxiter    = OPTIONS.maxiter;     end
if isfield(OPTIONS,'normb');         normb      = OPTIONS.normb;       end

gamma = pars.gamma;

if  (printyes)
    fprintf('\n *****************************************************');
    fprintf('******************************************');
    fprintf('\n \t   Semismooth Newton method for solving the subproblem');
    fprintf('\n ****************************************************');
    fprintf('*******************************************');
    fprintf('\n  iter cg_iter   lstep      grad_res         obj       cg_tol    time');
end


res_list = zeros(maxiter,1);

%% *************** the parameters for line search ***************

lsopts.c1 = 1e-4;

lsopts.c2 = 0.9;

lsopts.stepop = 1;

lsopts.printyes = 0;

cg_tol = 0.1;

cg_max = 10;   

cg_iter = 0; 

lstep = 1;

%% ********************** Main Loop *******************************

tstart = clock;

[fval,proxh,gradPsi] = fgrad(U,AtU,pars,Zk,Xk,AXk);

for iter=1:maxiter
    
    res = -gradPsi;

    norm_res = norm(res,'fro');
    
    res_list(iter) = norm_res;
    
    rnorm_res = norm_res/max(1,normb);
    
    %% ************* to check the stopping condition **************
    
    if (rnorm_res<=tol_SNCG)
     
        return;
        
    elseif((iter>=10)&& abs(norm_res-res_list(iter-5))/(1+abs(norm_res))<=1e-7)
   
        return;
        
    end
  
    ttime = etime(clock,tstart);
    
    if (printyes)
        
        fprintf('\n %3.0d     %2.0d     %3.2e     %3.2e    %5.6e      %3.2e    %3.2f',iter,cg_iter,lstep,rnorm_res,fval,cg_tol,ttime);
    end 
    
    cg_tol = min(cg_tol,norm_res^(1.2));
    
    ZkAtU = Zk + gamma*AtU;
    
 %% ************* calculate the generalized Newton direction *************
    
    tauj = min(1e-6,rnorm_res);
    
    [dir,cg_iter] = cg(res,res,pars,ZkAtU,Xk,cg_tol,cg_max,tauj);
    
    curve = sum(dot(gradPsi,dir));
    
    Atdir = Atmap(Xk,dir);
    
  %% ******************** Wolfe lsearch *************************  

   [U,AtU,fval,gradPsi,proxh,lstep] = Wolfe_cg(U,AtU,fval,curve,dir,Atdir,pars,Zk,Xk,AXk,lsopts,1e-6);
  
end
 




