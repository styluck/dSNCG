%% *****************************************************************
%  filename: Wolfe_cg
%
%% ****************************************************************

function [Unew,AtUnew,fnew,gnew,proxh,alp]= Wolfe_cg(Uk,AtUk,fk,curva0,dk,Atdk,pars,Zk,Xk,AXk,lsopts,tol)

stepop = lsopts.stepop;

printyes = lsopts.printyes;

maxit = ceil(log(1/tol)/log(2));

c1 = lsopts.c1;

c2 = lsopts.c2;

alpconst = 0.5;

%% ******************************************************************
for iter = 1:maxit
    
    if (iter==1)
        
        alp = 1; LB = 0; UB = 1;
    else
        
        alp = alpconst*(LB+UB);
    end
    
    Unew = Uk + alp*dk;
    
    AtUnew = AtUk + alp*Atdk;
    
    [fnew,proxh] = fgrad(Unew,AtUnew,pars,Zk,Xk,AXk);

    Vnew = proxh - Xk;
    
    curva = sum(dot(Vnew,Atdk));
    
    if (iter==1)
        
        gLB = curva0;  gUB = curva;
        
        if (sign(gLB)*sign(gUB)>0)

            if (printyes); fprintf('|'); end
                      
               gnew = -AXk + Amap(Xk,proxh);   
         
            return;
        end        
    end
    
    if (abs(curva)<c2*abs(curva0))&&(fnew-fk-c1*alp*curva0<1e-8/max(1,abs(fnew)))
        
        if ((stepop==1) || (stepop==2 && abs(curva)<tol))
            
            if (printyes); fprintf('|'); end
  
               gnew = -AXk + Amap(Xk,proxh);   
 
            return;
        end
    end
    
    if (sign(curva)*sign(gUB) < 0)
        
        LB = alp; gLB = curva;
        
    elseif (sign(curva)*sign(gLB)< 0)
        
        UB = alp; gUB = curva;
        
    end
    
    if (alp<1) && (printyes)
        fprintf('\n iter = %2d, ------line search value------------\n',iter);
        fprintf('\n ------alp = %2.2f, LQ = %4.3e, LQ0 = %4.3e',alp,fnew,fk);
    end
    
end
                
 gnew = -AXk + Amap(Xk,proxh);   

%% ****************************************************************