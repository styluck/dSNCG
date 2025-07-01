%% ****************************************************************
%  filename: dSNCG_BB
%% ****************************************************************

function [U,Xnew,fnew,gap_new,iter] = dSNCG_BB(U,Xk,gradk,fk,hfun,retr,pars,OPTIONS,K,n,r)

%%
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

%% ******************* parameters for SNDir ***********************

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
        
        [fnew,gnew] = objfun_QAP(Xnew,K, n,r, 1);
        
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
      
end  
    
%      Qp = Q; 
%      
%      Q = 0.85*Qp + 1; 
%     
%      Cval = (0.85*Qp*Cval + Fnew)/Q;    