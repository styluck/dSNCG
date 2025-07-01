%% **********************************************************************
%  filename: cg
%
%% **********************************************************************
%% conjugate gradient method for solving the system of linear equations
% 
%   H(X) = B  with  H = tauj + AVA'
%
%  where H is a linear operator from S^m to S^m
%
%  X0£ºthe starting point 
%
%  res£ºthe residual of LX = B at X0, i.e., res = B - L(X0) 
%%
%% **********************************************************************************  

function [x,iter,solve_ok] = cg(b,res,pars,U,Xk,tol,maxit,tauj)

m = size(b,1);

normb = norm(b,'fro');         % Here b is an m times m symmetric matrix

if ~exist('tol','var'); tol=1e-2*normb; end
 
if ~exist('maxit','var'); maxit = 10; end

resnrm = zeros(maxit+1,1);

solve_ok = 1;
 
tiny = 1.0e-16;
 
stagnate_check = 20;

%% ****************** Initialization part ***********************

d = zeros(m);   Hd = zeros(m);

x = zeros(m);   % such a starting point is crucial !!!

r = res;   err = norm(r,'fro');    
 
resnrm(1) = err;  

z = r;   rho_old = sum(sum(r.*z)); % 
 
tau_old = norm(z,'fro');
 
theta_old = 0;

%% ******************* Main Loop ********************************

for iter = 1:maxit
   
   [Hz,sigma] = Jacobi(z,pars,U,Xk,tauj);
        
   if (abs(sigma)<tiny)  %% now z=0 since H is positive definite
        
        x = res;
        
        solve_ok = -1;
        
        return;
        
    else
        
        alpha = rho_old/sigma;   %% this is the step-size
        
        r = r - alpha*Hz;
    end
    
    u = r; 
    
    normu = norm(r,'fro');
    
    theta = normu/tau_old;  
    
    c = 1/sqrt(1+theta^2);
            
    tau = c*normu;

    gamma = (c*theta_old)^2;
    
    eta = alpha/(1+theta^2);
    
    d = gamma*d + eta*z;    % the scaled direction
    
    Hd = gamma*Hd + eta*Hz;
    
    x = x + d;   
      
   %% ****************** stopping conditions  *******************
    
    res = res - Hd;
    
    err = norm(res,'fro');
    
    resnrm(iter+1) = err;

    if (err<tol)
        
        return; 
    end
    
    if (iter > stagnate_check) && (iter > 10)
        
        ratio = resnrm(iter-9:iter+1)./resnrm(iter-10:iter);
        
        if (min(ratio) > 0.997) && (max(ratio) < 1.003)
            
            solve_ok = -1;
  
            return;
        end
    end
    
  %% ************** update the new direction *********************
   
    rho = sum(sum(r.*u));
    
    beta = rho/rho_old;    % ||r||^2/||rold||^2
    
    z = u + beta*z;    
    
    rho_old = rho;
    
    tau_old = tau;
    
    theta_old = theta;
    
end
 
if (iter == maxit);  solve_ok = -2;  end
 
end    

%%%%%%%%%%%%%%%%%%%%%% End of conjugate_gradient.m  %%%%%%%%%%%%%%%%%%%%%%%      