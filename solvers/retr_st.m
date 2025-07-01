function Xnew = retr_st(Xk,Vk,r,method)

% The retraction function of Stiefel manifold
% supports for 4 types of retraction method, include:
% Exponential mapping (1),
% QR factorization (2),
% QR decomposition (3), and
% Cayley transform (4)
% Inputs:
%   Xk: the original point on Stiefel manifold
%   VK: the tangent vector in T_{Xk}M
%   method: the retraction method selected, default: Cayley
% Output:
%   Xnew: the new point on Stiefel manifold


if nargin < 3 || isempty(method)
    method = 4;
end

if method == 1
    I = eye(r);
    XkVk = Xk'*Vk;    
    [Q,~]=qr(-Vk+Xk*XkVk,0);    
    Z = zeros(r);
    Xnew = [Xk, Q] * expm([ -XkVk, -R'; R, Z])*[I; Z];
    
elseif method == 2
    I = eye(r);
    Xnew = (Xk + Vk)/(sqrt(I + Vk'*Vk));
    
elseif method == 3
    [Q,RR] = qr(Xk+Vk,0);
    diagRR = sign(diag(RR)); 
    ndr = diagRR < 0;
    if nnz(ndr)>0
       Q = Q*spdiags(diagRR,0,r,r);
    end
    Xnew = Q;
    
elseif method == 4
    
    n = size(Xk,1);
    
    I = eye(n);
    
    ImmXXtVXt = (I-0.5*Xk*Xk')*(Vk*Xk');
    
    W = ImmXXtVXt - ImmXXtVXt';
    
    Xnew =linsolve(I-0.5*W, Xk+0.5*Vk);
    
elseif method == 5
   
    PY = Xk+Vk;
    
    PYq = PY'*PY;
    
    PYq = 0.5*(PYq+PYq');
   
    [U,SIGMA,S] = eig(PYq);   
   
    SIGMA =diag(SIGMA);    
    
    Xnew = PY*(U*diag(sqrt(1./SIGMA))*S');

%       [U,SIGMA,S] = svd(PY);   
%     
%       SIGMA = diag(SIGMA);    
%      
%       Xnew = PY*(U*diag(sqrt(1./SIGMA))*S');
end
