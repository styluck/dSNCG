%% **************************************************************
%  filename: Moreau_h
%
%% ************ the envelope of the function ***********************
%
%  to compute the Moreau-envelop of the function 

%   h(X):=<E,max(0,-X)>   where E is a matrix of ones
%  
%   Zsol = P_{beta}h(Z):= argmin{(0.5/beta)*||X-Z||_F^2 + h(X)}

%   Mgam_h = (0.5/beta)*||Zsol-Z||_F^2 + h(Zsol)

%% ****************************************************************

function [Zsol,Mgam_h] = Moreau_h(Z,beta)

Zsol = min(Z+beta,max(Z,0));

if nargout>=2
    
    hval = -sum(Zsol(Zsol<0)); %sum(sum(max(0,-Zsol)));  
    
    Mgam_h = (0.5/beta)*norm(Zsol-Z,'fro')^2 + hval;
    
end