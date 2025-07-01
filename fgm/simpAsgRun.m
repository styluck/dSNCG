function wsRun = simpAsgRun(tagSrc, tagAlg, varargin)
% Run graph matching algorithm on the CMU Motion data set.
%
% Input
%   tagSrc  -  source type, 1 | 2 | 3
%   tagAlg  -  algorithm type, 1 | 2 | ...
%   varargin
%     save option
%
% Output
%   wsRun
%     prex  -  name
%     Me    -  mean, nAlg x nBin
%     Dev   -  standard deviation, nAlg x nBin
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 01-25-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 03-04-2013

% save option
prex = cellStr('cmum', 'tagSrc', tagSrc, 'tagAlg', tagAlg);
[svL, path] = psSv(varargin, ...
                   'prex', prex, ...
                   'subx', 'run', ...
                   'fold', 'cmum/asg/run');

% parameters for generating src
[tag, gaps, PFs, nIn] = cmumAsgPair(tagSrc);

% parameters for algorithms
[parAlgs, algs] = gmPar(tagAlg);

% dimension
nBin = length(gaps);
nReps = cellDim(PFs, 2);
nAlg = length(parAlgs);

% per gap
[Me, Dev, ObjMe, ObjDev, cpuTime] = zeross(nAlg, nBin);
prCIn('bin', nBin, 1);
for iBin = 1 : nBin
    prC(iBin);
    
    %src
    pFs = PFs{iBin}(:,1);
    wsSrc = cmumAsgSrc(tag, pFs, nIn);
    asgT = wsSrc.asgT;
    
    % feature
    parG = st('link', 'del');
    parF = st('smp', 'n', 'nBinT', 4, 'nBinR', 3);
    wsFeat = cmumAsgFeat(wsSrc, parG, parF);
    [gphs, XPs, Fs] = stFld(wsFeat, 'gphs', 'XPs', 'Fs');
    
    % affinity
    parKnl = st('alg', 'cmum');
    [KP, KQ] = conKnlGphPQU(gphs, parKnl);
    K = conKnlGphKU(KP, KQ, gphs);
    Ct = ones(size(KP));

    % undirected graph -> directed graph
    gphDs = gphU2Ds(gphs);
    KQD = [KQ, KQ; KQ, KQ];

    % per algorithm
    for iAlg = 1 : nAlg
        % parameter
        pars = parAlgs{iAlg};
        tstart = clock;
        if strcmpi(algs{iAlg}, 'fgm') || strcmpi(algs{iAlg}, 'fgm-u')
            asg = fgmU(KP, KQ, Ct, gphs, asgT, pars{:});
            asg.obj = asg.X(:)' * K * asg.X(:);
            
        elseif strcmpi(algs{iAlg}, 'fgm-d')
            asg = fgmD(KP, KQD, Ct, gphDs, asgT, pars{:});
            asg.obj = asg.X(:)' * K * asg.X(:);
            
        elseif strcmpi(algs{iAlg}, 'EPPGSN')
            asg = gm_dSNCG(K, asgT, pars{:});
            
        elseif strcmpi(algs{iAlg}, 'SEPPG')
            asg = gm_Mpenalty(K, asgT, pars{:});
            
        elseif strcmpi(algs{iAlg}, 'Spenalty-BB')
            asg = gm_Spenalty(K, asgT, pars{:});
            
        elseif strcmpi(algs{iAlg}, 'ALM')
            asg = gm_QAPORTHALEQQ(K, asgT, pars{:});
            
        elseif strcmpi(algs{iAlg}, 'ALMre')
            asg = gm_QAPORTHALEQQ0(K, asgT, pars{:});
            
        elseif strcmpi(algs{iAlg}, 'ep4orth')
            asg = gm_ep4orth(K, asgT, pars{:});
            
        else
            asg = gm(K, Ct, asgT, pars{:});
        end
        
        % objective
        Time(iAlg, 1) = etime(clock,tstart);
        Acc(iAlg, 1) = asg.acc;
        Obj(iAlg, 1) = asg.obj;
    end
    
    % mean & deviation
    Me(:, iBin) = mean(Acc, 2);
    Dev(:, iBin) = std(Acc, 0, 2);
    ObjMe(:, iBin) = mean(Obj, 2);
    ObjDev(:, iBin) = std(Obj, 0, 2);
    cpuTime(:, iBin) = mean(Time,2);
end
prCOut(nBin + 1);

% store
wsRun.prex = prex;
wsRun.Me = Me;
wsRun.Dev = Dev;
wsRun.ObjMe = ObjMe;
wsRun.ObjDev = ObjDev;
wsRun.cpuTime = cpuTime;
% save
if svL > 0
    save(path, 'wsRun');
end

prOut;
