%This script creates the data for the phase retrieval instance: This is a 
%psd matrix with dimension n = 25,921. It has exact rank 250, but its 
%effective rank R = 5. The matrix is an approximate solution to a phase 
%retrieval SDP. The solution is obtained by the approach described in 
%[YHC2015] and using the data from [YUTC2017]. See our referneces for 
%more details.
%
%This code generates a random matrix with same size and singular values as 
%the dataset. The original dataset is not shared due to COPYRIGHT.
%
%	[TYUC2019] J.A. Tropp, A. Yurtsever, M. Udell and V. Cevher.  
%	Streaming Low-Rank Matrix Approximation with an Application to
%	Scientific Simulation. 
%
%	[YUTC2017] A. Yurtsever, M. Udell, J.A. Tropp and V. Cevher. Sketchy 
%	Decisions: Convex Low-Rank Matrix Optimization with Optimal Storage. 
%
%   [YHC2015] A. Yurtsever, Y.P. Hsieh and V. Cevher. Scalable Convex 
%   Methods for Phase Retrieval. 
%
%	Coded by: Alp Yurtsever
%	Ecole Polytechnique Federale de Lausanne, Switzerland.
%   Laboratory for Information and Inference Systems, LIONS.
%   contact: alp.yurtsever@epfl.ch
%   Last modified: July 16, 2018 (Modified from Nys-SKETCH)
%
% SKETCHv1.1
% Copyright (C) 2018 Laboratory for Information and Inference Systems
% (LIONS), Ecole Polytechnique Federale de Lausanne, Switzerland.
% This code is a part of SKETCH toolbox. 
% Please read COPYRIGHT before using this file.

%% Fix the seed
clearvars
rng(0,'twister');

%% Test setup
load data/PhaseRetrievalSingVals
info.U = randn(info.m, size(info.Z,1)) + 1i*randn(info.m, size(info.Z,1));
info.U = orth(info.U);
save data/PhaseRetrieval