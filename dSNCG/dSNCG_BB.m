%% ****************************************************************
%  filename: dSNCG_BB
%% ****************************************************************

function [U,Xnew,fnew,gap_new,iter] = dSNCG_BB(U,Xk,gradk,fk,hfun,retr,pars,OPTIONS,A,B,r)
% conjugate gradient method for solving the system of linear equations
% 
%   H(X) = B  with  H = tauj + AVA'
%
%  where H is a linear operator from S^m to S^m
%
%  X0 -- the starting point 
%
%  res -- the residual of LX = B at X0, i.e., res = B - L(X0) 

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
 

%
if isfield(OPTIONS,'maxiter');     maxiter    = OPTIONS.maxiter;     end
if isfield(OPTIONS,'printyes');    printyes   = OPTIONS.printyes;    end
if isfield(OPTIONS,'Vtol');        Vtol       = OPTIONS.Vtol;        end
if isfield(OPTIONS,'Ftol');        Ftol       = OPTIONS.Ftol;        end


min_gam = pars.min_gam; 

max_gam = pars.max_gam;  

eta = 0.1;

rhoLS = 1.0e-4;

rho = pars.rho;

gamma = pars.gamma;

AXk = 2*eye(r);   % Amap(Xk,Xk);

% ******************* parameters for SNDir ***********************

OPTIONS_SNDir.printyes = 0;

OPTIONS_SNDir.maxiter = 10;  %% this is the best !!

OPTIONS_SNDir.normb = 2*sqrt(r);

OPTIONS_SNDir.tol = 1e-8;

%% ****************************************************************

gapk = hfun(Xk);

AtU = Atmap(Xk,U);

Ck = ProjXZ(Xk,gradk); 

Zk = Xk - gamma*Ck;

Fk_list = zeros(maxiter,1);

%% ********************* main loop ********************************

%Q = 1;  

Cval = fk + rho*gapk; 

for iter = 1:maxiter
    
    Uold = U; AtU_old = AtU;
   
    nls = 1;
    
    while 1   %% to search a desired step-size
                        
        [U,AtU,Ysol,itSNCG] = SNewton(U,AtU,pars,OPTIONS_SNDir,Zk,Xk,AXk);
 
        Vk = Ysol-Xk;
        
        normVk = norm(Vk,'fro');
        
        sqnormVk = normVk^2;        
 
        Xnew = retr(Xk,Vk);
        
        [fnew,gnew] = objfun_QAP(Xnew,A,B);
        
        gap_new = hfun(Xnew);
        
        Fnew = fnew + rho*gap_new;
        
        if Fnew<=Cval-(0.5/gamma)*rhoLS*sqnormVk || nls>=5
     
            break;
        end
        
        gamma = eta*gamma;
        
        pars.gamma = gamma;
        
        nls = nls + 1;
        
        U = Uold;  AtU = AtU_old;
    end
         
    Fk_list(iter) = Fnew;
            
    if normVk <=Vtol
        
        return;
        
    elseif (normVk<=5*Vtol && iter>=10 && abs(Fnew-Fk_list(iter-9))/(1+abs(Fnew))<=Ftol)

        return;
    end
        
    Cnew = ProjXZ(Xnew,gnew);      %% do not negelect !!
    
 %% *************** to estimate the step-size via BB *****************
    
    DetaX = Xnew - Xk;  
     
    DetaY = Cnew - Ck; 
           
    DetaXY = abs(sum(dot(DetaX,DetaY)));
    
    gamma1 = norm(DetaX,'fro')^2/DetaXY;
    
    gamma2 = DetaXY/norm(DetaY,'fro')^2;
    
    gamma = max(min(min(gamma1,gamma2),max_gam),min_gam);
    
    pars.gamma = gamma;
    
    if (printyes)
        
        fprintf('\n %3d    %3d       %3.2e      %5.4e     %3.2e    %3.2e    %2.1e   %2d',iter,itSNCG,sqnormVk,fnew,gamma,rho,gapk,nls);
        
    end
       
    Xk = Xnew;  Ck = Cnew;
            
    Zk = Xk - gamma*Ck; 
    
    gapk = gap_new;  
    
    if iter<=5
        
       Cval = max(Fk_list(1:iter));
    else
            
       Cval = max(Fk_list(iter-5+1:iter));
    end
    
    Vtol = max(0.99*Vtol,1.0e-3);    
end  
    
%     Qp = Q; 
%      
%      Q = 0.85*Qp + 1; 
%     
%      Cval = (0.85*Qp*Cval + Fnew)/Q;    
%
%    temp_vec = Fk_list(iter-5+1:iter);
        
%    IdxK = find(temp_vec<Fnew+1000/(iter^(2.1)));
        
%    Cval = max(temp_vec(IdxK));