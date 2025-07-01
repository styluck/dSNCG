function [pars, algs] = gmPar(tag, lmxitr, objX)
% Obtain parameters for graph matching algorithm.
%
% Input
%   tag     -  type of pair, 1 | 2
%                1 : ga, sm, smac, ipfp1, ipfp2, rrwm, fgmU
%                2 : ga, sm, smac, ipfp1, ipfp2, rrwm, fgmU, fgmD
%
% Output
%   pars    -  parameters for each algorithm, 1 x nAlg (cell), 1 x nPari (cell)
%   algs    -  algorithm name, 1 x nAlg (cell)
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011
%   modify  -  Feng Zhou (zhfe99@gmail.com), 05-06-2013
%   modify  -  Q Tian & S Pan & L Xiao (lxiao@scut.edu.cn), 08-20-2021
if nargin < 3 || isempty(objX)
    objX = 1; % quadratic (1) or quartic (2)
end

if nargin < 2 || isempty(lmxitr)
    lmxitr = 10;
end

[pars, algs] = cellss(1, 100);

% graph matching with asymmetric edge feature

nAlg = 0;
if tag == 2
    %-----------------
    % IPFP-U
    nAlg = nAlg + 1;
    parIni = st('alg', 'unif', 'nor', 'doub');
    parPosC = st('alg', 'none');
    parPosD = st('alg', 'ipfp', 'deb', 'n');
    pars{nAlg} = {parIni, parPosC, parPosD};
    algs{nAlg} = 'IPFP-U';
    
    % IPFP-S
    nAlg = nAlg + 1;
    parIni = st('alg', 'sm', 'top', 'eigs');
    parPosC = st('alg', 'none');
    parPosD = st('alg', 'ipfp', 'deb', 'n');
    pars{nAlg} = {parIni, parPosC, parPosD};
    algs{nAlg} = 'IPFP-S';
    
    % RRWM
    nAlg = nAlg + 1;
    parIni = st('alg', 'unif', 'nor', 'unit');
    parPosC = st('alg', 'rrwm');
    parPosD = st('alg', 'hun');
    pars{nAlg} = {parIni, parPosC, parPosD};
    algs{nAlg} = 'RRWM';
    
    % FGM-U
    nAlg = nAlg + 1;
    parFgmS = st('nItMa', 100, 'nAlp', 101, 'thAlp', 0, 'deb', 'n', 'ip', 'n');
    pars{nAlg} = {parFgmS};
    algs{nAlg} = 'FGM-U';
    
    % FGM-D
    nAlg = nAlg + 1;
    parFgmA = st('nItMa', 100, 'nAlp', 101, 'deb', 'n', 'ip', 'n', 'lamQ', .5);
    pars{nAlg} = {parFgmA};
    algs{nAlg} = 'FGM-D';
    
    %-----------------
    % dSNCG_QAPBB
   nAlg = nAlg + 1;
    pardSNCG = st('maxiter', 2000, 'printyes', 0, 'gtol', 1e-3, 'ftol', 1e-8,...
        'lmxitr',lmxitr,'objX', objX);
    pars{nAlg} = {pardSNCG};
    algs{nAlg} = 'dSNCG-QAPBB';
    
    % Mpenalty_BB
    nAlg = nAlg + 1;
    parMpenalty = st('maxiter', 2000, 'printyes', 0, 'gtol', 1e-3, 'objtol', 1e-8,...
        'lmxitr',lmxitr,'objX', objX);
    pars{nAlg} = {parMpenalty};
    algs{nAlg} = 'SEPPG';
    
    % Spenalty_BB
    nAlg = nAlg + 1;
    % parSpenalty = st('maxiter', 2000, 'printyes', 0, 'gtol', 1e-3, 'ftol', 1e-8);
    pars{nAlg} = {parMpenalty};
    algs{nAlg} = 'Spenalty-BB';
    
    % QAPORTHALEQQ_QUAD
    nAlg = nAlg + 1;
    parQAPORTHALEQQ = st('record', 0, 'tol', 1e-3, 'tolsub', 1e-3, 'omxitr', 100, ...
        'mxitr', 100, 'xtol', 1e-5, 'gtol', 1e-5, 'ftol', 1e-8, ...
        'tau', 1e-3, 'objX', objX,'lmxitr',lmxitr);
    pars{nAlg} = {parQAPORTHALEQQ};
    algs{nAlg} = 'ALM';
    
    % QAPORTHALEQQ_QUADre
    nAlg = nAlg + 1;
    parQAPORTHALEQQ0 = st('record', 0, 'tol', 1e-3, 'tolsub', 1e-3, 'omxitr', 100, ...
        'mxitr', 100, 'xtol', 1e-5, 'gtol', 1e-5, 'ftol', 1e-8, ...
        'tau', 1e-3, 'objX', objX,'lmxitr',lmxitr);
    pars{nAlg} = {parQAPORTHALEQQ0};
    algs{nAlg} = 'ALMre';
    % ep4orth
    nAlg = nAlg + 1;
    parep4orth = st('sigma', 1e-3, 'omega2', 1.05, 'omega3', 0.98, 'altmethod', inf,...
        'lmxitr',lmxitr,'objX', objX);
    pars{nAlg} = {parep4orth};
    algs{nAlg} = 'ep4orth';
    
elseif tag == 1
    % dSNCG_QAPBB with objX = 1
    nAlg = nAlg + 1;
    pardSNCG = st('maxiter', 2000, 'printyes', 0, 'gtol', 1e-3, 'ftol', 1e-8,...
        'lmxitr',lmxitr,'objX', 2);
    pars{nAlg} = {pardSNCG};
    algs{nAlg} = 'EPPGSN';
    
elseif tag == 3
    % FGM-D
    nAlg = nAlg + 1;
    parFgmA = st('nItMa', 100, 'nAlp', 101, 'deb', 'n', 'ip', 'n', 'lamQ', .5);
    pars{nAlg} = {parFgmA};
    algs{nAlg} = 'FGM-D';
    
    % dSNCG_QAPBB
    nAlg = nAlg + 1;
    pardSNCG = st('maxiter', 2000, 'printyes', 0, 'gtol', 1e-3, 'ftol', 1e-8,...
        'lmxitr',lmxitr,'objX', 2);
    pars{nAlg} = {pardSNCG};
    algs{nAlg} = 'EPPGSN';
    
    % Mpenalty_BB
    nAlg = nAlg + 1;
    parMpenalty = st('maxiter', 2000, 'printyes', 0, 'gtol', 1e-3, 'objtol', 1e-8,...
        'lmxitr',lmxitr,'objX', objX);
    pars{nAlg} = {parMpenalty};
    algs{nAlg} = 'SEPPG';
%     
%     % Spenalty_BB
%     nAlg = nAlg + 1;
%     % parSpenalty = st('maxiter', 2000, 'printyes', 0, 'gtol', 1e-3, 'ftol', 1e-8);
%     pars{nAlg} = {parMpenalty};
%     algs{nAlg} = 'Spenalty-BB';
%     
    % QAPORTHALEQQ_QUAD
    nAlg = nAlg + 1;
    parQAPORTHALEQQ = st('record', 0, 'tol', 1e-3, 'tolsub', 1e-3, 'omxitr', 100, ...
        'mxitr', 100, 'xtol', 1e-5, 'gtol', 1e-5, 'ftol', 1e-8, ...
        'tau', 1e-3, 'objX', objX,'lmxitr',lmxitr);
    pars{nAlg} = {parQAPORTHALEQQ};
    algs{nAlg} = 'ALM';
    
    % QAPORTHALEQQ_QUADre
    nAlg = nAlg + 1;
    parQAPORTHALEQQ0 = st('record', 0, 'tol', 1e-3, 'tolsub', 1e-3, 'omxitr', 100, ...
        'mxitr', 100, 'xtol', 1e-5, 'gtol', 1e-5, 'ftol', 1e-8, ...
        'tau', 1e-3, 'objX', objX,'lmxitr',lmxitr);
    pars{nAlg} = {parQAPORTHALEQQ0};
    algs{nAlg} = 'ALMre';
    
    % ep4orth
    nAlg = nAlg + 1;
    parep4orth = st('sigma', 1e-3, 'omega2', 1.05, 'omega3', 0.98, 'altmethod', inf,...
        'lmxitr',lmxitr,'objX', objX);
    pars{nAlg} = {parep4orth};
    algs{nAlg} = 'ep4orth';
    
elseif tag == 4
    % ep4orth
    nAlg = nAlg + 1;
    parep4orth = st('sigma', 1e-3, 'omega2', 1.05, 'omega3', 0.98, 'altmethod', inf,...
        'lmxitr',lmxitr,'objX', objX);
    pars{nAlg} = {parep4orth};
    algs{nAlg} = 'ep4orth';
end

% outs
pars(nAlg + 1 : end) = [];
algs(nAlg + 1 : end) = [];

% [EOF]