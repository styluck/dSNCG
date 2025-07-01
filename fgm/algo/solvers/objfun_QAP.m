%% ***************************************************************
%  filename: objfun_QAP
%
%% *****************************************************************
%  Compute the augmented Lagrangian and its gradient  
%
%% ***************************************************************
%% where Ck,Zk and AXk are constant matrices !!
%%

function [fval,gradf]= objfun_QAP(X,K, n, r, mtd)
if nargin < 5
    mtd = 1;
end

if nargin < 3
    [n, r] = size(X);
end
X = X(:);
if mtd == 1
    KX = K*X;
    fval = - X'*KX; % objective func is maximizing
    if nargout >= 2
    %     [n1, r1] = size(KX);
    %     fprintf('n=%2.2f, r=%2.2f, n1=%2.2f, r1=%2.2f\n',n,r,n1,r1);
        % grad = AXB + A*Xq*B';

        gradf = reshape( - (KX + K'*X), n, r);

    end
elseif mtd == 2
    Xq = X.^2;
    KX = K*Xq;

    fval = - Xq'*KX; % objective func is maximizing

    if nargout >= 2
    %     [n1, r1] = size(KX);
    %     fprintf('n=%2.2f, r=%2.2f, n1=%2.2f, r1=%2.2f\n',n,r,n1,r1);
        % grad = AXB + A*Xq*B';

        gradf = reshape( -2 * (KX + K'*Xq).*X, n, r);

    end
end
