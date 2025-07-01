function prIni()
% propmter for random initialization
global Iter
if isempty(Iter)
    Iter = 1;
end
randstate = 100+(Iter-1)*9867;
randn('state',double(randstate));
rand('state',double(randstate));
state = rng;
Iter = Iter + 1;
% [EOF]