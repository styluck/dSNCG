%% ***************************************************************
%  filename: Jacob
%
%% ****************************************************************
%  Calculate the direction: tauj*dir + A*V*A'dir
%
%  V:   the Clarke Jacobian of A Prox_{rho*gamma}h(Zk+gamma*AtU)
%  dir: an r times r matrix
%% **************************************************************

function [Hd,sigma,Atd] = Jacobi(d,pars,ZkAtU,Xk,tauj)

%tauj = 1.0e-7;

gamma = pars.gamma;

rho = pars.rho;

Atd = Atmap(Xk,d);      % Here d is a r times r matrix

%W = zeros(n,r);        % Clarke Jacobian of Proxh

W=(ZkAtU>=0 &ZkAtU<=-gamma*rho);

%% ********************** Jacobian of Proxh ***********************

WAtd = W.*Atd;

AWAtd = Amap(Xk,WAtd);

Hd = tauj*d + gamma*AWAtd;

sigma = tauj*norm(d,'fro')^2 + gamma*sum(dot(WAtd,WAtd));

end