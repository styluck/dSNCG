%% ***************************************************************
%  filename: MEpen_fgrad
%
%% *****************************************************************
%  Compute the Moreau-envelop of nonsmooth exact penalty  
%
%% ***************************************************************
%% where Ck,Zk and AXk are constant matrices !!
%%

function [fval,grad]= MEpen_fgrad(X,A,B,rho)

lmbda = 1e+0;

Xq = X.^2;

%% ********************** objective value and gradient ***************

AXB = A'*Xq*B;  

[proxh,Mh] = Moreau_h(X,lmbda);

fval = sum(dot(Xq,AXB)) + rho*Mh;

if nargout>=2
    
   grad = 2*(AXB+(A*Xq*B')).*X + (rho/lmbda)*(X-proxh);
  
%   grad = 2*(AXB+(A*Xq*B')).*X -(1+p)*rho*max(0,-X).^p;
end