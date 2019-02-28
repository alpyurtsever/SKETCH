function [sd] = Test_single_LowRankPlusNoise(n,r,R,noiseLevel,field,model)
%This script runs an experimental instance using a matrix from the low rank 
%plus noise class. See our reference paper for more details.
%	Inputs: n       -> dimension of matrix (n x n)
%           r		-> target rank
%			R		-> effective rank
%			noiseLevel	->  'LowNoise' (low noise) 
%							'MedNoise' (medium noise) 
%							'HiNoise'  (high noise)
% 			field 	-> field for the test matrices: 'complex' or 'real'.
%           model   -> model for the test matrices: 'Gauss','Sparse','SSRFT'
%	Output: sd		-> structure with following fields
%						* ErrBest_S2 ,ErrBest_S1, ErrBest_Sinf
%						* ID, twosketch, threesketch
%						* info (struct with more info about size and storage)
%						* twosketch (struct with info of error for [Upa2016] & [TYUC2017])
%						* threesketch (struct with info of error for [TYUC2019])
%						* hmtsketch (struct with info of error for [HMT2011])
%						* wlrtsketch (struct with info of error for [WLRT2008])
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
%   Last modified: June 22, 2018
%
% SKETCHv1.1
% Copyright (C) 2018 Laboratory for Information and Inference Systems
% (LIONS), Ecole Polytechnique Federale de Lausanne, Switzerland.
% This code is a part of SKETCH toolbox. 
% Please read COPYRIGHT before using this file.

dataname = ['data/LowRank',noiseLevel,'_n',num2str(n),'_R',num2str(R),'.mat'];

if ~exist(dataname,'file')
    error(['Data file not found: Before go to the data folder and run ' ... 
    '''CreateData.m'' file to generate a dataset for noisy low-rank' ...
    'model before running this script.']);
end
load(dataname)

ErrBest_S2 = norm(singVals(r+1:end));
ErrBest_S1 = norm(singVals(r+1:end),1);
ErrBest_Sinf = norm(singVals(r+1:end),inf);

%% Run the generic test file
Test_generic_Type1;

end