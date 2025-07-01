%% ***************************************************************
%  filename: Epen_fgrad
%
%% *****************************************************************
%  Compute the augmented Lagrangian and its gradient
%
%% ***************************************************************
%% where Ck,Zk and AXk are constant matrices !!
%%

function [fval,grad,objf]= Epen_fgrad(X,A, rho, n, r, mtd)

if nargin < 6
    mtd = 1;
end

X = X(:);

Mh = sum(X(X<0).^2); %sum(sum(max(0,-X).^2));

if mtd == 1
    AXB = A*X;
    objf = - sum(dot(X, AXB));
    fval = objf + rho*Mh;
    if nargout>=2
        
        grad = reshape( - (AXB + A'*X) + (2*rho)*min(0,X), n, r);
        
        %   grad = 2*(AXB+(A*Xq*B')).*X -(1+p)*rho*max(0,-X).^p;
    end
    
elseif mtd == 2
    Xq = X.^2;
    %% ********************** objective value and gradient ***************
    
    AXB = A*Xq;
    objf = - sum(dot(Xq, AXB));
    fval = objf + rho*Mh;
    
    if nargout>=2
        
        grad = reshape( - 2*(AXB + A'*Xq).*X + (2*rho)*min(0,X), n, r);
        
        %   grad = 2*(AXB+(A*Xq*B')).*X -(1+p)*rho*max(0,-X).^p;
    end
end