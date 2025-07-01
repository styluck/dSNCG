%% ***************************************************************
%  filename: MEpen_fgrad
%
%% *****************************************************************
%  Compute the Moreau-envelop of nonsmooth exact penalty
%
%% ***************************************************************
%% where Ck,Zk and AXk are constant matrices !!
%%

function [fval,grad,objf]= MEpen_fgrad(X,A,rho,n,r,mtd)

if nargin < 6
    mtd = 1;
end

lambda = 0.05;

X = X(:);

[proxh,Mh] = Moreau_h(X,lambda);

if mtd == 1
    
    AXB = A*X;
    
    objf = -sum(dot(X,AXB)); 
    
    fval = objf + rho*Mh;
    
    if nargout>=2        
        grad = reshape(-(AXB+A'*X)+ (rho/lambda)*(X-proxh), n, r);
    end
    
elseif mtd == 2
    
    Xq = X.^2;
    
   %% ********************** objective value and gradient ***************
    
    AXB = A*Xq;
    
    objf = -Xq'*AXB;
    
    fval = objf + rho*Mh;
    
    if nargout>=2
        
        grad = reshape(-2*(AXB+A'*Xq).*X+(rho/lambda)*(X-proxh),n,r);
    end
end

% [EOF]