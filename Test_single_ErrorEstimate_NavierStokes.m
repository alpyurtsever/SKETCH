function [sd] = Test_single_ErrorEstimate_NavierStokes(field, model)
%This file implements the experiments for the scatter plots in Fig.6 of the
%main reference [TYUC2019], which presents the sampling distributions of
%the approximation error and the error estimator. 
%
%We use various values of k (we fix s=2k+1) and various values of error 
%sketch size q (3,5,7,10 are used in the tests, only 5 and 10 are used in 
%plots). 
%
%   See our reference paper for the detailed explanation of the
%   sketching procedure and the arithmetic, communication and storage
%   costs.
%
%	[TYUC2019] J.A. Tropp, A. Yurtsever, M. Udell and V. Cevher.  
%	Streaming Low-Rank Matrix Approximation with an Application to
%	Scientific Simulation. 
%
%   Coded by: Alp Yurtsever
%   Ecole Polytechnique Federale de Lausanne, Switzerland.
%   Laboratory for Information and Inference Systems, LIONS.
%   contact: alp.yurtsever@epfl.ch
%   Last modified: January 30, 2019
%
% SKETCHv1.1
% Copyright (C) 2018 Laboratory for Information and Inference Systems
% (LIONS), Ecole Polytechnique Federale de Lausanne, Switzerland.
% This code is a part of SKETCH toolbox. 
% Please read COPYRIGHT before using this file.

%% Generate test setup and/or load data

load('data/NavierStokes-StreamVelocity.mat');

%% Run the generic test file

if strcmp(field,'real')
    alpha = 1;
    beta = 1;
elseif strcmp(field,'complex')
    alpha = 0;
    beta = 2;
else
    error('''field'' should be ''real'' or ''complex''.')
end

[m,n] = size(A);
if m < n
    A = A';
    [m,n] = size(A);
end

kSweep = [8, 16, 32, 48, 96, 128];
lKs = length(kSweep);
qSweep = [3,5,7,10];
sSweep = nan(lKs,1);
spectrumLowRankApprox = zeros(lKs,lKs);
spectrumFixedRankApprox = zeros(lKs,lKs,lKs);
spectrumExact = singVals;
ErrLowRank = nan(lKs,1);
ErrFixedRank = nan(lKs,lKs);

ErrLowRankEstimate = {};
ErrFixedRankEstimate = {};
Theta = {};
S = {};
for q = qSweep
    Theta{end+1} = Gauss(q,m,field);                %#ok
    S{end+1} = Theta{end}*A;                        %#ok
    ErrLowRankEstimate{end+1} = nan(lKs,1);         %#ok
    ErrFixedRankEstimate{end+1} = nan(lKs,lKs);     %#ok
end

for k = kSweep
    s = 2*k + 1;
    sSweep(k) = s;
    
    %% Create the Sketch
    myThreeSketch = ThreeSketch(model, m, n, k, s, field);
    myThreeSketch.LinearUpdate(A);
    
    %% Low Rank Approximation
    [Q,W,P] = myThreeSketch.LowRankApprox();
    spectrumLowRankApprox(k,1:k) = svd(W);
    
    Ahat = Q*W*P';
    ErrLR = norm(A - Ahat,'fro');
    ErrLowRank(k) = ErrLR;
    
    %% Error Sketch
    for t = 1:length(qSweep)
        q = qSweep(t);
        Sapprox = ((Theta{t}*Q)*W)*P';
        ErrLREst = (1/sqrt(beta*q)) * norm(S{t} - Sapprox,'fro');
        ErrLowRankEstimate{t}(k) = ErrLREst;        %#ok
    end
    
    %% Fixed Rank Approximation
    for r = 1:k
        [U,Sigma,V] = myThreeSketch.FixedRankApprox(r);
        Ahat_r = U*Sigma*V';
        spectrumFixedRankApprox(k,r,1:r) = diag(Sigma);
        
        ErrFR = norm(A - Ahat_r,'fro');
        ErrFixedRank(k,r) = ErrFR;
        
        for t = 1:length(qSweep)
            q = qSweep(t);
            Sapprox = ((Theta{t}*U)*Sigma)*V';
            ErrFREst = (1/sqrt(beta*q)) * norm(S{t} - Sapprox,'fro');
            ErrFixedRankEstimate{t}(k,r) = ErrFREst;    %#ok
        end
    end
    
    %% Print intermediate results
    % fprintf('k = %d, s = %d \n', k, s);
end

sd.kSweep = kSweep;
sd.qSweep = qSweep;
sd.sSweep = sSweep;
sd.spectrumApprox = spectrumLowRankApprox;
sd.spectrumFixedRankApprox = spectrumFixedRankApprox;
sd.spectrumExact = spectrumExact;
sd.ErrLowRank = ErrLowRank;
sd.ErrLowRankEstimate = ErrLowRankEstimate;
sd.ErrFixedRank = ErrFixedRank;
sd.ErrFixedRankEstimate = ErrFixedRankEstimate;
sd.ID = [datestr(now,30),num2str(randi([100,999]))];

% savepath = ['results/',mfilename,'/'];
% if ~exist(savepath,'dir'), mkdir(savepath); end
% filename = datestr(now,30);
% save([savepath,filename],'-v7.3','sd');

end
