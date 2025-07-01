%% ***************************************************************
%  filename: objfun_QAP
%
%% *****************************************************************
%  Compute the augmented Lagrangian and its gradient  
%
%% ***************************************************************
%% where Ck,Zk and AXk are constant matrices !!
%%

function [fval,gradf]= objfun_QAP(X,A,B)

Xq = X.^2;

AXB = A'*Xq*B;       

fval = sum(dot(Xq,AXB));

if nargout>=2
  
  % gradf = AXB + A*Xq*B';
    
   gradf = 2*(AXB+(A*Xq*B')).*X;
end