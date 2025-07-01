%% ****************************************************************
%  filename: PG_BB
%% ****************************************************************

function [Xnew,iter,objf] = PG_BB(Xk,retr,pars,OPTIONS,A)

%%
if isfield(OPTIONS,'maxiter');     maxiter    = OPTIONS.maxiter;     end
if isfield(OPTIONS,'printyes');    printyes   = OPTIONS.printyes;    end
if isfield(OPTIONS,'Vtol');        Vtol       = OPTIONS.Vtol;        end
if isfield(OPTIONS,'Ftol');        Ftol       = OPTIONS.Ftol;        end

min_gam = pars.min_gam;

max_gam = pars.max_gam;

eta = 0.3;

rhoLS = 1.0e-4;

rho = pars.rho;

gamma = pars.gamma;

%% ****************************************************************

[Fold,gradF] = Epen_fgrad(Xk,A,rho);

grad = ProjXZ(Xk,gradF);

Fk_list = zeros(maxiter,1);

%% ********************* main loop ********************************

%Q = 1;  

Cval = Fold;

for iter = 1:maxiter
    
    nls = 1;
    
    while 1   %% to search a desired step-size
        
        Vk = -gamma*grad; 
        
        normVk = norm(Vk,'fro');
        
        sqnormVk = normVk^2;
        
        Xnew = retr(Xk,Vk);
        
        [Fnew,gradF_new,objf] = Epen_fgrad(Xnew,A,rho);
        
        if (Fnew<= Cval-(0.5/gamma)*rhoLS*sqnormVk) %||(nls>=5)
            
            break;
        end
        
        gamma = eta*gamma;
        
        nls = nls + 1;
    end

    Fk_list(iter) = Fnew;
    
    if normVk<=Vtol
        
        return;
        
    elseif (normVk<=5*Vtol && iter>=10 && abs(Fnew-Fk_list(iter-9))/(1+abs(Fnew))<=Ftol)
        
        return;
    end
    
    grad_new = ProjXZ(Xnew,gradF_new);
    
%% *************** to estimate the step-size via BB *****************
    
    DetaX = Xnew - Xk;  DetaY = grad_new - grad;
    
    DetaXY = abs(sum(dot(DetaX,DetaY)));
    
    gamma1 = norm(DetaX,'fro')^2/DetaXY;
    
    gamma2 = DetaXY/norm(DetaY,'fro')^2;
    
    gamma = max(min(min(gamma1,gamma2),max_gam),min_gam);
    
    if (printyes)&&mod(iter,1)==0
        
        fprintf('\n %3d     %3.2e      %5.4e     %3.2e    %3.2e    %2d',iter,sqnormVk,Fnew,gamma,rho,nls);
        
    end
    
    Xk = Xnew;  grad = grad_new;

    if iter<=5        
        Cval = max(Fk_list(1:iter));
    else
        Cval = max(Fk_list(iter-5+1:iter));
    end
    
    Vtol = max(0.99*Vtol,1.0e-3);
end


    
    %     Qp = Q;
    %
    %     Q = 0.85*Qp + 1;
    %
    %     Cval =(0.85*Qp*Cval + Fnew)/Q;
