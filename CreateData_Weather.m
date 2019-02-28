%This script creates the data for the Weathe instance: The real (m x n) 
%matrix MinTemp contains the minimum temperature recorded at each of 
%m = 19264 stations on each of n = 7305 days. The first 10 singular 
%values decay by two orders of magnitude, while the rest of the spectrum
%has medium polynomial decay. This is typical for measured data. See 
%our reference paper for more details.
%
%This code generates a random matrix with same size and singular values as 
%the dataset. The original dataset is not shared due to COPYRIGHT.
%
%	[TYUC2019] J.A. Tropp, A. Yurtsever, M. Udell and V. Cevher.  
%	Streaming Low-Rank Matrix Approximation with an Application to
%	Scientific Simulation. 
%
%   Coded by: Alp Yurtsever
%   Ecole Polytechnique Federale de Lausanne, Switzerland.
%   Laboratory for Information and Inference Systems, LIONS.
%   contact: alp.yurtsever@epfl.ch
%   Last modified: July 16, 2018
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
load data/WeatherSingVals
m = 19264;
n = 7305;
[U,~] = qr(randn(m,n),0);
[V,~] = qr(randn(n,n));
A = U*(spdiags(singVals,0,n,n)*V');
clearvars U V m n;
save data/Weather.mat