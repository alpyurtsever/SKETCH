%This script batches the experimental instances into the cluster.
%Experiments are carried out in parallel. 
%
%	[TYUC2019] J.A. Tropp, A. Yurtsever, M. Udell and V. Cevher.  
%	Streaming Low-Rank Matrix Approximation with an Application to
%	Scientific Simulation. 
%
%   Coded by: Alp Yurtsever
%   Ecole Polytechnique Federale de Lausanne, Switzerland.
%   Laboratory for Information and Inference Systems, LIONS.
%   contact: alp.yurtsever@epfl.ch
%   Last modified: June 22, 2018
%
% SKETCHv1.1
% Copyright (C) 2018 Laboratory for Information and Inference Systems
% (LIONS), Ecole Polytechnique Federale de Lausanne, Switzerland.
% This code is a part of SKETCH toolbox. 
% Please read COPYRIGHT before using this file.

%% Beginning of the test

MonteCarloNum = 20;
model = 'Sparse';

c = parcluster;

ClusterInfo.setMemUsage('16gb');

experiment          = {};

experiment{end+1}   = @() Test_parfor_Weather(10,'real',model,MonteCarloNum);
experiment{end+1}   = @() Test_parfor_NavierStokes(10,'StreamVelocity','real',model,MonteCarloNum);
experiment{end+1}   = @() Test_parfor_MaxCut(1,'real',model,MonteCarloNum);
experiment{end+1}   = @() Test_parfor_MaxCut(14,'real',model,MonteCarloNum);
experiment{end+1}   = @() Test_parfor_PhaseRetrieval(1,'complex',model,MonteCarloNum);
experiment{end+1}   = @() Test_parfor_PhaseRetrieval(5,'complex',model,MonteCarloNum);

%% 
j = {};
for tt = 1:length(experiment)
    fprintf('Batching experiment No:%d...',tt);
    c.batch(experiment{tt},1,{},'pool',MonteCarloNum,...
        'CurrentFolder','YOUR-WORKING-DIRECTORY'); 	% You need to Set this folder as your working directory on your cluster
    fprintf(' Started..\n');
end
