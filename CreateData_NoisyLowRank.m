%This script draws a random input matrix instances for Low-Rank + Noise
%model. See our reference paper for more details.
%
%	[TYUC2019] J.A. Tropp, A. Yurtsever, M. Udell and V. Cevher.  
%	Streaming Low-Rank Matrix Approximation with an Application to
%	Scientific Simulation. 
%
%	Coded by: Alp Yurtsever
%	Ecole Polytechnique Federale de Lausanne, Switzerland.
%	Laboratory for Information and Inference Systems, LIONS.
%	contact: alp.yurtsever@epfl.ch
%	Last modified: July 16, 2018 (Mod. from PRACTICALSKETCH/Nys-SKETCH)
%
% SKETCHv1.1
% Copyright (C) 2018 Laboratory for Information and Inference Systems
% (LIONS), Ecole Polytechnique Federale de Lausanne, Switzerland.
% This code is a part of SKETCH toolbox. 
% Please read COPYRIGHT before using this file.

%% Fix the seed
rng(0,'twister');

%% Test setup
clearvars

xi = [1e-4, 1e-2, 1e-1];
sname = {'LowNoise','MedNoise','HiNoise'};

n = 1000;

Rsweep = [5,10,20];

for R = Rsweep
    for t = 1:length(xi)
        
        % Generate the matrix A
        A = zeros(n,n);
        A(1:R,1:R) = diag(ones(R,1));
        G = (randn(n,n) + 1i*randn(n,n));
        W = G*G';
        clearvars G;
        A = A + xi(t)/n.*W;
        clearvars W;
        
        %% Find the ground truth (best rank-r approximation)
        [~, Sigmar, ~] = svd(A);
        
        singVals = diag(Sigmar);
        clearvars Sigmar;

        %% Save data sparse
        dataname = ['data/LowRank',sname{t},'_n',num2str(n),'_R',num2str(R),'.mat'];
        save(dataname)
                
    end
end