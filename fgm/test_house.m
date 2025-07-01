% Testing the performance of different graph matching methods on the CMU House dataset.
% This is the same function used for reporting (Fig. 4) the first experiment (Sec 5.1) in the CVPR 2013 paper.
%
% History
%   create   -  Feng Zhou (zhfe99@gmail.com), 09-01-2011
%   modify   -  Feng Zhou (zhfe99@gmail.com), 05-08-2013
%   modify   -  Y Qian & S Pan & L Xiao (lxiao@scut.edu.cn), 08-20-2021

clear variables;

restoredefaultpath;
addPath;

prSet(3);

%% save flag
svL = 2; % change svL = 1 if you want to re-run the experiments.
tagAlg = 3;
objX = 1;
lmxitr = 30;
%% algorithm parameter
[~, algs] = gmPar(tagAlg,30,1); % load params

%% run 1 (perfect graphs, no noise)
clc
tagSrc = 1;
[~, val1s] = cmumAsgPair(tagSrc); %load target source
wsRun1 = cmumAsgRun(tagSrc, tagAlg, 'svL', svL);

%% run 2 (randomly remove nodes)
tagSrc = 2;
[~, val2s] = cmumAsgPair(tagSrc);
wsRun2 = cmumAsgRun(tagSrc, tagAlg, 'svL', svL);

%% cal run data
[~, algs] = gmPar(tagAlg); % load params
algs(6) = {'SNCG'};
algs(7) = {'SEPPG'};
algs(8) = {'EPPGSN'};
[Me1, Dev1, ObjMe1, ObjDev1] = stFld(wsRun1, 'Me', 'Dev', 'ObjMe', 'ObjDev');
om1 = mean(ObjMe1);
oma1 = max(max(ObjMe1));
ObjMe1 = ObjMe1./oma1;
ObjDev1 = ObjDev1./oma1;

[Me2, Dev2, ObjMe2, ObjDev2] = stFld(wsRun2, 'Me', 'Dev', 'ObjMe', 'ObjDev');
om2 = mean(ObjMe2);
oma2 = max(max(ObjMe2));
ObjMe2 = ObjMe2./oma2;
ObjDev2 = ObjDev2./oma2;

%% show accuracy & objective
rows = 1; cols = 3;

Ax = iniAx(1, rows, cols, [500 * rows, 500 * cols]);

shCur(Me1, Dev1, 'ax', Ax{1}, 'dev', 'y');
xticks = 1 : 3 : size(Me1, 2);
setAxTick('x', '%d', xticks, val1s(xticks));
set(gca, 'ylim', [.7 1.05], 'ytick', .1 : .1 : 1, 'xlim', [.5, size(Me1, 2) + .5]);
grid on
axis square;
xlabel('baseline');
ylabel('accuracy');

shCur(ObjMe1, ObjDev1, 'ax', Ax{2}, 'dev', 'y');
xticks = 1 : 3 : size(ObjMe1, 2);
setAxTick('x', '%d', xticks, val1s(xticks));
set(gca, 'ylim', [.7 1.05], 'ytick', .1 : .1 : 1, 'xlim', [.5, size(Me1, 2) + .5]);
grid on
axis square;
xlabel('baseline');
ylabel('objective');

shCur(Me1, Dev1, 'ax', Ax{3}, 'dev', 'n', 'algs', algs);
% set(Ax{3}, 'visible', 'off');

%%
rows = 1; cols = 3;

Ax = iniAx(1, rows, cols, [500 * rows, 500 * cols]);

shCur(Me2, Dev2, 'ax', Ax{1}, 'dev', 'y');
xticks = 1 : 3 : size(Me2, 2);
setAxTick('x', '%d', xticks, val1s(xticks));
set(gca, 'ylim', [.7 1.05], 'ytick', .1 : .1 : 1, 'xlim', [.5, size(Me1, 2) + .5]);
grid on
axis square;
xlabel('baseline');
ylabel('accuracy');

shCur(ObjMe2, ObjDev2, 'ax', Ax{2}, 'dev', 'y');
xticks = 1 : 3 : size(ObjMe1, 2);
setAxTick('x', '%d', xticks, val1s(xticks));
set(gca, 'ylim', [.7 1.05], 'ytick', .1 : .1 : 1, 'xlim', [.5, size(Me1, 2) + .5]);
grid on
axis square;
xlabel('baseline');
ylabel('objective');

shCur(Me1, Dev1, 'ax', Ax{3}, 'dev', 'n', 'algs', algs);
% set(Ax{5}, 'visible', 'off');
% cla;
