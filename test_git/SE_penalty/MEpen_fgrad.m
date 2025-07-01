%% ***************************************************************
%  filename: MEpen_fgrad
%
%% *****************************************************************
%  Compute the Moreau-envelop of nonsmooth exact penalty
%
%% ***************************************************************
%% where Ck,Zk and AXk are constant matrices !!
%%

function [fval,grad]= MEpen_fgrad(X,A, rho)

lambda = 1e+0;


%% ********************** objective value and gradient ***************

[proxh,Mh] = Moreau_h(X,lambda);

[fobj,gobj] = objfun(X,A);

fval = fobj + rho*Mh;

if nargout>=2
    
    grad = gobj + (rho/lambda)*(X-proxh);
    
end