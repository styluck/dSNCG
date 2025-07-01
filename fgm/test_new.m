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


%% save flag
svL = 2; % change svL = 1 if you want to re-run the experiments.
tagAlg = 3;
objX = 1;
lmxitr = 20;
prSet(tagAlg);
%% algorithm parameter
[~, algs] = gmPar(tagAlg,lmxitr,objX); % load params

%% run 1 (house, perfect graphs, no noise)
tagSrc = 1;
[~, val1s] = cmumAsgPair(tagSrc); %load target source
wsRun1 = simpAsgRun(tagSrc, tagAlg, 'svL', svL);

%% run 2 (house, randomly remove 5 nodes)
tagSrc = 2;
[~, val2s] = cmumAsgPair(tagSrc);
wsRun2 = simpAsgRun(tagSrc, tagAlg, 'svL', svL);
% %% run 3 (hotel, perfect graphs, no noise)
% tagSrc = 3;
% [~, val3s] = cmumAsgPair(tagSrc); %load target source
% wsRun3 = simpAsgRun(tagSrc, tagAlg, 'svL', svL);
% 
% %% run 4 (hotel, randomly remove 5 nodes)
% tagSrc = 4;
% [~, val4s] = cmumAsgPair(tagSrc);
% wsRun4 = simpAsgRun(tagSrc, tagAlg, 'svL', svL);
% 
% if exist('wsRun1','var') && exist('wsRun2','var')
%     save(['wsRun',datestr(now,'mmdd'),'.mat'],'wsRun1','wsRun2');
% end

%% plotings
scrsz = get(0,'ScreenSize');
[mks, cls] = genMkCl;
len = size(val1s, 2);

for i = 1:4
    if exist(sprintf('wsRun%d',i),'var')
        eval(sprintf('wsRun = wsRun%d;',i));
        eval(sprintf('vals = val%ds;',i));
        MaxObj = max(max(wsRun.ObjMe));
        wsRun.ObjMe = wsRun.ObjMe/MaxObj;
        
        h = figure('Position',[100 scrsz(4)*.22 scrsz(3)*.8 scrsz(4)*.4]);
        subplot(1,3,1)
        hold on
        for j = 1:size(wsRun.Me,1)
            plot(wsRun.Me(j,:)','Marker',mks{j}, 'Color', cls{j}, 'LineWidth', 1)
        end 
        
        xticks = 1 : 3 : len;
        setAxTick('x', '%d', xticks, vals(xticks));
        set(gca, 'ylim', [.7 1.05], 'ytick', .1 : .1 : 1, 'xlim', [.5, len + .5]);
        grid on
        axis square;
        xlabel('Sequence gap');
        ylabel('Accuracy');
        legend(algs,'location','southwest')
        
        subplot(1,3,2)
        hold on
        for j = 1:size(wsRun.ObjMe,1)
            plot(wsRun.ObjMe(j,:)','Marker',mks{j}, 'Color', cls{j}, 'LineWidth', 1)
        end 
        
        xticks = 1 : 3 : len;
        setAxTick('x', '%d', xticks, vals(xticks));
        set(gca, 'ylim', [.7 1.05], 'ytick', .1 : .1 : 1, 'xlim', [.5, len + .5]);
        grid on
        axis square;
        xlabel('Sequence gap');
        ylabel('Objective');
        legend(algs,'location','southwest')
        
        subplot(1,3,3)
        hold on
        for j = 1:size(wsRun.cpuTime,1)
            plot(wsRun.cpuTime(j,:)','Marker',mks{j}, 'Color', cls{j}, 'LineWidth', 1)
        end 
        
        xticks = 1 : 3 : len;
        setAxTick('x', '%d', xticks, vals(xticks));
        set(gca, 'ylim', [0 200], 'ytick', 0 : 20 : 200, 'xlim', [.5, len + .5]);
        grid on
        axis square;
        xlabel('Sequence gap');
        ylabel('CPU time');
        legend(algs,'location','northwest')
    end
end


% [EOF]