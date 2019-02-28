%This script creates the data for the Navier Stokes Stream Velocity data. 
%This data is courtesy of Beverley McKeon and Sean Symon. 
%The real (m x n) matrix StreamVel contains streamwise velocities at 
%m = 10738 points for each of n = 5001 time instants. The first 20 singular 
%values of the matrix decay by two orders of magnitude, and the rest of 
%the spectrum exhibits slow exponential decay. This is typical for physical 
%models. See our reference paper for more details.
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
load data/StreamVelSingVals
m = 10738;
n = 5001;
[U,~] = qr(randn(m,n),0);
[V,~] = qr(randn(n,n));
A = U*(spdiags(singVals,0,n,n)*V');
clearvars U V m n;
save data/NavierStokes-StreamVelocity.mat