%% ****************************************************************
%  filename: Proj_orth
%
%% ****************************************************************
%% to compute argmin ||Z-G||_F^2 s.t. Z^TZ=I
%%

function X = Proj_orth(G)
    
   [U,D,V] = svd(G,'econ');   % U is n times r and V is r times r
   
   X = U*V';

end