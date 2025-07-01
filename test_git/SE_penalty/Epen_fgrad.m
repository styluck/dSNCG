%% ***************************************************************
%  filename: Epen_fgrad
%
%% *****************************************************************
%  Compute the augmented Lagrangian and its gradient
%
%% ***************************************************************
%% where Ck,Zk and AXk are constant matrices !!
%%

function [fval,grad,fobj]= Epen_fgrad(X,A, rho)



%% ********************** objective value and gradient ***************

Mh = sum(X(X<0).^2); %sum(sum(max(0,-X).^2));

[fobj,gobj] = objfun(X,A);

fval = fobj + rho*Mh;

if nargout>=2
    
    grad = gobj - (2*rho)*min(0,X);
end