%This script runs an experimental instance using a real valued Navier-Stokes  
%dataset. The results are used to plot Fig.4 and Fig.5 in [TYUC2019]. 
%See our reference paper for more details.
%
%   [WLRT2008] F. Woolfe, E. Liberty, V. Rokhlin and M. Tygert.
%   A Fast Randomized Algorithm for the Approximation of Matrices.
%
%   [HMT2011] N. Halko, P.G. Martinsson and J.A. Tropp. Finding Structure
%   with Randomness: Probabilistic algorithms for constructing approximate 
%   matrix decompositions,
%
%   [Upa2016] J. Upadhyay, Fast and Space-Optimal Low-Rank Factorization in 
%   the Streaming Model with Application in Differential Privacy.
%
%	[TYUC2017] J.A. Tropp, A. Yurtsever, M. Udell and V. Cevher.  
%	Practical Sketching Algorithms for Low-Rank Matrix Approximation. 
%
% 	[TYUC2019] J.A. Tropp, A. Yurtsever, M. Udell and V. Cevher.  
% 	Streaming Low-Rank Matrix Approximation with an Application to
% 	Scientific Simulation. 
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

rng(0,'twister');

load('data/NavierStokes-StreamVelocity.mat');

field = 'real';

%% Beginning of the test

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

kSweep = 1:128;
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
ErrZeroEstimate = {};
Theta = {};
S = {};
for q = qSweep
    Theta{end+1} = Gauss(q,m,field);                                %#ok
    S{end+1} = Theta{end}*A;                                        %#ok
    ErrLowRankEstimate{end+1} = nan(lKs,1);                         %#ok
    ErrFixedRankEstimate{end+1} = nan(lKs,lKs);                     %#ok
    ErrZeroEstimate{end+1} = (1/sqrt(beta*q))*norm(S{end},'fro');   %#ok
end

for k = kSweep
    s = 2*k + 1;
    sSweep(k) = s;
    
    %% Create the Sketch
    myThreeSketch = ThreeSketch('Sparse', m, n, k, s, field);
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
    fprintf('k = %d, s = %d \n', k, s);
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
sd.ErrZeroEstimate = ErrZeroEstimate;
sd.ID = [datestr(now,30),num2str(randi([100,999]))];

savepath = ['results/',mfilename,'/'];
if ~exist(savepath,'dir'), mkdir(savepath); end
filename = datestr(now,30);
save([savepath,filename],'-v7.3','sd');
