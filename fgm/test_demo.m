% A demo comparison of different graph matching methods on the on CMU House dataset.
%
% Remark
%   The edge is directed and the edge feature is asymmetric.
% clc;
% clear;
% close all;

% addPath;
% make;

%% src params
tag = 'house';
pFs = [1 31]; % frame index
nIn = [30 30] - 5; % randomly remove 5 nodes
parKnl = st('alg', 'cmum'); % type of affinity: only edge distance

% algorithm parameter1
[pars, algs] = gmPar(2,10,2);

% src
wsSrc = cmumAsgSrc(tag, pFs, nIn, 'svL', 1);
asgT = wsSrc.asgT;

% feature
parG = st('link', 'del'); % Delaunay triangulation for computing the graphs
parF = st('smp', 'n', 'nBinT', 4, 'nBinR', 3); % not used, ignore it
wsFeat = cmumAsgFeat(wsSrc, parG, parF, 'svL', 1);
[gphs, XPs, Fs] = stFld(wsFeat, 'gphs', 'XPs', 'Fs');

% affinity
[KP, KQ] = conKnlGphPQU(gphs, parKnl);
K = conKnlGphKU(KP, KQ, gphs);
Ct = ones(size(KP));
scaler = 1; % 1; % 1e+04;
% pop out
fprintf('\n 1: dSNCG_QAPBB\n 2: Mpenalty_BB\n 3: Spenalty_BB\n 4: Wen & Yin\n');
fprintf(' 5~9: other solvers for graph matching\n');
fprintf(' 10:  epen4orth\n');
%% solvers
Slct = 2;% input('Select a solver(1~9):\n'); %8;
tic
switch Slct
    case 1  % dSNCG_QAPBB % 0.483
        Kppa = K*scaler;
        asg = gm_dSNCG(Kppa, asgT, pars{6}{:}); % PPASNCG
%         asg = gm_dSNCG(K, asgT, pars{6}{:}); % PPASNCG
        fprintf('best acc=%2.3f\n',asg.acc);
        fprintf('avg acc = %2.3f\n', mean(asg.accsq));
    case 2 % Mpenalty_BB % 0.692
        Kppa = K*scaler;
        asg = gm_Mpenalty(Kppa, asgT, pars{7}{:}); % PPASNCG
        fprintf('best acc=%2.3f\n',asg.acc);
        fprintf('avg acc = %2.3f\n', mean(asg.accsq));
    case 3 % Spenalty_BB % 0.741
        Kppa = K*scaler;
        asg = gm_Spenalty(Kppa, asgT, pars{8}{:}); % Wen & Yin StiefelGBB
        fprintf('best acc=%2.3f\n',asg.acc);
        fprintf('avg acc = %2.3f\n', mean(asg.accsq));
    case 4 % Wen & Yin StiefelGBB %0.447
        Kppa = K*scaler;
        asg = gm_QAPORTHALEQQ(Kppa, asgT, pars{9}{:}); % Wen & Yin StiefelGBB
        fprintf('best acc=%2.3f\n',asg.acc);   
        fprintf('avg acc = %2.3f\n', mean(asg.accsq));
    case 5 % IPFP-U
        asg = gm(K, Ct, asgT, pars{1}{:});
        fprintf('best acc=%2.3f\n', asg.acc);
        
    case 6 % IPFP-S
        asg = gm(K, Ct, asgT, pars{2}{:});
        fprintf('best acc=%2.3f\n', asg.acc);
    case 7 % RRWM
        asg = gm(K, Ct, asgT, pars{3}{:});
        fprintf('best acc=%2.3f\n', asg.acc);
    case 8 % FGM-U
        asg = fgmU(KP, KQ, Ct, gphs, asgT, pars{4}{:});
        fprintf('best acc=%2.3f\n', asg.acc);
    case 9 % FGM-D
        % undirected graph -> directed graph (for FGM-D)
        gphDs = gphU2Ds(gphs);
        KQD = [KQ, KQ; KQ, KQ];
        asg = fgmD(KP, KQD, Ct, gphDs, asgT, pars{5}{:});
        fprintf('best acc=%2.3f\n', asg.acc);
    case 10 
        Kppa = K*scaler;
        asg = gm_ep4orth(Kppa, asgT, pars{7}{:}); % ep4orth
        fprintf('best acc=%2.3f\n',asg.acc);
        fprintf('avg acc = %2.3f\n', mean(asg.accsq));
end
toc
%% show correspondence result
rows = 1; cols = 1;

Ax = iniAx(1, rows, cols, [400 * rows, 900 * cols], 'hGap', .1, 'wGap', .1);
parCor = st('cor', 'ln', 'mkSiz', 7, 'cls', {'y', 'b', 'g'});
shAsgImg(Fs, gphs, asg , asgT, parCor, 'ax', Ax{1}, 'ord', 'n');
% title(sprintf('result of %s. acc: %2.2f',asg.alg,asg.acc));
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,'fig1_1','-dpdf','-r0')
% [EOF]