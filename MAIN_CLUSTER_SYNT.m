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

n = 1000;
r = 10;
field = 'complex';
MonteCarloNum = 20;

c = parcluster;

ClusterInfo.setMemUsage('4gb');

experiment          = {};

for R = [5,10,20]
for model = {'Gauss','SSRFT','Sparse'}
% Low-Rank + noise
experiment{end+1}   = @() Test_parfor_LowRankPlusNoise(n,r,R,'HiNoise',field,model{1},MonteCarloNum);   %#ok
experiment{end+1}   = @() Test_parfor_LowRankPlusNoise(n,r,R,'MedNoise',field,model{1},MonteCarloNum);  %#ok
experiment{end+1}   = @() Test_parfor_LowRankPlusNoise(n,r,R,'LowNoise',field,model{1},MonteCarloNum);  %#ok
% % Low-Rank + poly decay
experiment{end+1}   = @() Test_parfor_PolyDecay(n,r,R,0.5,field,model{1},MonteCarloNum);                %#ok
experiment{end+1}   = @() Test_parfor_PolyDecay(n,r,R,1,field,model{1},MonteCarloNum);                  %#ok
experiment{end+1}   = @() Test_parfor_PolyDecay(n,r,R,2,field,model{1},MonteCarloNum);                  %#ok
% % Low-Rank + exp decay
experiment{end+1}   = @() Test_parfor_ExpDecay(n,r,R,0.01,field,model{1},MonteCarloNum);                %#ok
experiment{end+1}   = @() Test_parfor_ExpDecay(n,r,R,0.1,field,model{1},MonteCarloNum);                 %#ok
experiment{end+1}   = @() Test_parfor_ExpDecay(n,r,R,0.5,field,model{1},MonteCarloNum);                 %#ok
end
end


%% 
j = {};
for tt = 1:length(experiment)
    fprintf('Batching experiment No:%d...',tt);
    c.batch(experiment{tt},1,{},'pool',MonteCarloNum,...
        'CurrentFolder','YOUR-WORKING-DIRECTORY'); 	% You need to Set this folder as your working directory on your cluster
    fprintf(' Started..\n');
end
