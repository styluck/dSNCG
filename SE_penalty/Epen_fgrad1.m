xsdj%% ***************************************************************
%  filename: Epen_fgrad
%
%% *****************************************************************
%  Compute the augmented Lagrangian and its gradient  
%
%% ***************************************************************
%% where Ck,Zk and AXk are constant matrices !!
%%

function [fval,grad,objf]= Epen_fgrad(X,A,B,rho)


Xq = X.^2;

%% ********************** objective value and gradient ***************

AXB = A'*Xq*B;  

Mh = sum(X(X<0).^2); %sum(sum(max(0,-X).^2));

objf = sum(dot(Xq,AXB));

fval = objf + rho*Mh;

if nargout>=2
    
   grad = 2*(AXB+(A*Xq*B')).*X +(2*rho)*min(0,X);

end