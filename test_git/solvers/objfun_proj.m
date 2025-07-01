
function [fval,gradf]= objfun_proj(X,A)
% objective function for proj

fval = sum(dot(X,X))-2*sum(dot(X,A));

if nargout>=2
    gradf = 2*(X - A);
end

% [EOF]