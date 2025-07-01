%% ***************************************************************
%  filename: fgrad
%
%% *****************************************************************
%  Compute the function and gradient values of Psi  
%
%  Psi(U)=gamma||U||_F^2-rho*M_{gamma*rho}vtheta(Zk+gamma*AtU) 
%
%% ***************************************************************
%% where Ck,Zk and AXk are constant matrices !!
%%

function [fval,proxh,gradPsi]= fgrad(U,AtU,pars,Zk,Xk,AXk)

gamma = pars.gamma;

rho = pars.rho;

%% **************** The value of Psi and its gradient **************

tempU = Zk + gamma*AtU;

[proxh,Mh] = Moreau_h(tempU,gamma*rho);

%fval = (0.5*gamma)*norm(AtU-Ck,'fro')^2 - rho*Mh;

fval = 2*gamma*norm(U,'fro')^2 - rho*Mh + (0.5/gamma)*norm(Xk-Zk,'fro')^2; %%

if nargout>=3

   gradPsi = -AXk + Amap(Xk,proxh);
   
end



