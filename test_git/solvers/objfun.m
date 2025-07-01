function [fval,gradf]= objfun(X,A)

slct = 1;
switch slct
    
    case 1 % objective function for proj
        fval = sum(dot(X,X))-2*sum(dot(X,A));
        
        if nargout>=2
            gradf = 2*(X - A);
        end
        
    case 2 % objective function for onmf
        C = (X'*X)\X'*A;
        C = max(C,0);
        C_C = 2*(C*C');
        A_C =  2*(A*C');
        
        fval = -sum(dot(A_C,X))+sum(dot(C,(X'*X)*C));
        
        if nargout>=2
            gradf = (X*C_C-A_C);
        end
end
% [EOF]