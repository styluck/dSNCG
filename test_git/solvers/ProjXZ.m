%% ***************************************************************
%   filename: ProjXZ
%% ***************************************************************
% used to compute the gradient under Euclidean metric 

function Proj_mat = ProjXZ(X,Z)

XtZ = X'*Z;

tempXZ = 0.5*(XtZ + XtZ');

Proj_mat = Z - X*tempXZ; 