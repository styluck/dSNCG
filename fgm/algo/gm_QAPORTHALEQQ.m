function asg = gm_QAPORTHALEQQ(K, asgT, opts)
% Graph matching using Wen-Yin method 
%
% Input
%   K         -  global kernel matrix, nn x nn
%   asgT      -  ground-truth assignment (can be [])
% Output
%   asg       -  assignment
%     alg     -  'fgmU'
%     X       -  correspondence matrix, n1 x n2
%     acc     -  accuracy
%     obj     -  objective
% 
% Ref: Wen, Z., & Yin, W. (2013). 
% A feasible method for optimization with orthogonality constraints. 
% Mathematical Programming, 142(1), 397-434.

% function parameter
Ksnorm = svds(K,1);
lmxitr = opts.lmxitr; % 101

ha = tic;

[n, r] = size(asgT.X);

% params
acc = -100;
obj= -100;

% subopts.record = 0;
% subopts.mxitr = 100;
% subopts.xtol = 1e-5;
% subopts.gtol = 1e-5;
% subopts.ftol = 1e-8;

accsq = zeros(lmxitr, 1);
objsq = zeros(lmxitr, 1);

% algorithm: Wen-Yin solver
for i =1:lmxitr
    randstate = 100+(i-1)*9867;
    randn('state',double(randstate));
    rand('state',double(randstate));
    state = rng;
    
    X0 = randn(n);    
    X0 = orth(X0);
    %X0 = OptStiefelGBB(X0, @objfun_QAP,subopts,K,n,r,2);
    [X_tmp, out] = QAPORTHALEQQ_QUAD(X0, K, 0.1, opts, n,r);
    
    X_tmp = round(X_tmp);
    accsq(i) = matchAsg(X_tmp, asgT);
    objsq(i) = X_tmp(:)'*K*X_tmp(:);
% compare with ground-truth: the accuracy
    if objsq(i) > obj
        X = X_tmp;
        acc = accsq(i);
        obj = objsq(i);
    end
end


% store
asg.acc = acc;
asg.alg = 'ALM';
asg.X = X;
asg.obj = obj;
asg.tim = toc(ha);
asg.accsq = accsq;
asg.objsq = objsq;
end

% [EOF]