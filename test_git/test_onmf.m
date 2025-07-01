% function test_onmf
%--------------------------------------------------------------------------
% Orthogonal nonnegative matrix factorization (ONMF)
%
%                  min     \|A-XY^T\|_F^2
%                  s.t.    X'X=I, X>=0, Y>=0.
%
% Data available at http://www.cad.zju.edu.cn/home/dengcai/Data/data.html
%
%--------------------------------------------------------------------------
% Reference:
% B. Jiang, X. Meng, Z. Wen and X. Chen
% An Exact Penalty Approach For Optimization With Nonnegative Orthogonality
% Constraints
%
% Author: X. Meng, B. Jiang
% modified: Y. Qian, S. Pan, L. Xiao
% Version 1.0 .... 2021/1

%--------------------------------------------------------------------------

% add path
clc; clear; close all;
restoredefaultpath;
src = pwd;
addpath(genpath(src));

% choose examples
Problist = 1:9;

filesrc = strcat(pwd,filesep,'results');
if ~exist(filesrc, 'dir');     mkdir(filesrc);   end
filepath = strcat(filesrc, filesep, 'onmf');
if ~exist(filepath, 'dir');    mkdir(filepath);  end
strnum = '%6s %2d %+5.4e %+5.4e %+5.4e %+5.4e\n';
fprintf('\t\t\t\t  pur\t\t  ent\t\t time\t\tfeasi\n');

% load funcs
% sqhfun = @(x)sum(x(x<0).^2);      % sum(sum(max(0,-x).^2))
% hfun = @(x)abs(sum(x(x<0)));      % sum(sum(max(0,-x)))
% set params for algos
opts.maxiter = 1000;
opts.printyes = 0;
opts.gtol = 1.0e-3;
opts.ftol = 1.0e-8; % for dSNCG
opts.objtol = 1.0e-8; % for Spenalty/Mpenalty

for dsolver = 1
    
    switch dsolver
        case 1
            fprintf('------ result from dSNCG ------\n');
        case 2 
            fprintf('------ result from Spenalty_BB ------\n');
        case 3
            fprintf('------ result from Mpenalty_BB ------\n');
    end
    
    for dprob = Problist
        % get matrix A and ground truth
        switch dprob
            case 1
                data = load('2k2k.mat');
                A = data.fea;
                test_label = 10;
                true_ans = data.gnd+1;
            case 2
                data = load('Yale_32x32.mat');
                A = data.fea;
                test_label = 15;
                true_ans = data.gnd;
            case 3
                tra = load('TDT2-l10.mat');
                A = tra.A;
                true_ans = tra.true_ans;
                test_label = 10;
            case 4
                tra = load('TDT2-l20.mat');
                A = tra.A;
                true_ans = tra.true_ans;
                test_label = 20;
            case 5
                tra = load('TDT2-t10.mat');
                A = tra.A;
                true_ans = tra.true_ans;
                test_label = 10;
            case 6
                tra = load('TDT2-t20.mat');
                A = tra.A;
                true_ans = tra.true_ans;
                test_label = 20;
            case 7
                tra = load('Reu-t10.mat');
                A = tra.A;
                true_ans = tra.true_ans;
                test_label = 10;
            case 8
                tra = load('Reu-t20.mat');
                A = tra.A;
                true_ans = tra.true_ans;
                test_label = 20;
            case 9
                tra = load('News-t5.mat');
                A = tra.A;
                true_ans = tra.true_ans;
                test_label = 5;
        end
        
        % preprocessing
        [n_,m_] = size(A);
        sumA = full(sum(A));
        A2 = zeros(n_,m_);
        couu = 1;
        for i=1:m_
            if(sumA(i)~=0)
                A2(:,couu) = A(:,i);
                couu = couu+1;
            end
        end
        A = A2(:,1:couu-1);
        [n, m] = size(A);
        
        retr = @(x,v) retr_st(x,v,test_label,3);
        % get the initial point
        ptimetic = tic;
        [X0, ~, ~] = svds(A,test_label);
        %     for j=1:test_label
        %         neg = X0(:,j)<0;
        %         if norm(X0(neg,j),'fro') > norm(X0(~neg,j),'fro')
        %             X0(:,j) = -X0(:,j);
        %         end
        %     end
        %     X0 = max(X0,0);
        %     X0_norm = sqrt(sum(X0.*X0));
        %     X0 = X0./X0_norm;
        pretime = toc(ptimetic);
        
        % set parameters for algos
        switch dsolver
            case 1
                ctimetic = tic;
                [fobj1, gap, X] = dSNCG_QAPBB(X0,opts,A);
                ctime = toc(ctimetic);
            case 3
                ctimetic = tic;
                [fobj2,Xsol2] = Mpenalty_BB(X0,opts,A);
                ctime = toc(ctimetic);
            case 2
                ctimetic = tic;
                [fobj2,Xsol2] = Spenalty_BB(X0,opts,A);
                ctime = toc(ctimetic);
        end
        
        % record results, check purity and other criteria
        time = ctime + pretime;
        [~, wh] = max(X,[],2);
        Y = zeros(n,test_label);
        for i=1:n
            Y(i,wh(i)) = 1;
        end
        [pur, ent] = comp_pe(Y,true_ans,1e-8);
        feasi = norm(X'*X-eye(test_label),'fro');
        fprintf(strnum, 'case',dprob,pur,ent,time,feasi);
        
        
    end
end


% [EOF]