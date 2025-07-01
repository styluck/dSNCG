%% ***************************************************************
%   filename: CProjXZ
%% ***************************************************************
% used to compute the gradient under canonical metric 

function Proj_mat = CProjXZ(X,Z)

ZtX = Z'*X;

Proj_mat = Z - X*ZtX; 