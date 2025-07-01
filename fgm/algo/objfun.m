%% ***************************************************************
%  filename: ALM_fun
%
%% *****************************************************************
%  Compute the augmented Lagrangian and its gradient  
%
%% ***************************************************************
%% where Ck,Zk and AXk are constant matrices !!
%%

function [fval,gradf]= objfun(X,K, n, r)

if nargin < 3
    [n, r] = size(X);
end
vecX = X(:);

KX = K*vecX;

fval = - vecX'*KX; % objective func is maximizing

if nargout >= 2
%     [n1, r1] = size(KX);
%     fprintf('n=%2.2f, r=%2.2f, n1=%2.2f, r1=%2.2f\n',n,r,n1,r1);
    % grad = AXB + A*Xq*B';
    
    gradf = reshape( - (KX + K'*vecX) , n, r);
    
end
% ------------------------------------------------
% ------------------------------------------------

% [n, r] = size(X);
% 
% vecX = X(:);
% 
% Xq = vecX.^2;
% 
% KX = K*Xq;
% 
% fval = - sum(sum(Xq.*KX)); % objective func is maximizing
% 
% if nargout >= 2
%     
%     % grad = AXB + A*Xq*B';
%     
%     gradf = reshape(- 2*(KX + K'*Xq).*vecX, n, r);
%     
% end
end