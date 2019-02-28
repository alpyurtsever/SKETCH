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

c = parcluster;

ClusterInfo.setMemUsage('16gb');

experiment          = {};

experiment{end+1} = @() Test_parfor_ErrorEstimate_NavierStokes('real','Sparse',1000);
experiment{end+1} = @() Test_ErrorEstimate_Figure_4_5('real','Sparse',1000);

%% 
j = {};
for tt = 1:length(experiment)
    fprintf('Batching experiment No:%d...',tt);
    c.batch(experiment{tt},1,{},'pool',100,...
        'CurrentFolder','YOUR-WORKING-DIRECTORY'); 	% You need to Set this folder as your working directory on your cluster
    fprintf(' Started..\n');
end
