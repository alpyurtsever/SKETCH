%This script creates the data for the Max-Cut instance: This is a 
%real-valued psd matrix with dimension n = 2000, and its effective rank
%R = 14. The matrix is an approximate solution to the MAXCUT SDP for the 
%sparse graph G40. See our reference paper for more details.
%
%	[TYUC2019] J.A. Tropp, A. Yurtsever, M. Udell and V. Cevher.  
%	Streaming Low-Rank Matrix Approximation with an Application to
%	Scientific Simulation. 
%
%	Coded by: Alp Yurtsever
%	Ecole Polytechnique Federale de Lausanne, Switzerland.
%	Laboratory for Information and Inference Systems, LIONS.
%	contact: alp.yurtsever@epfl.ch
%	Last modified: July 16, 2018 (Modified from Nys-SKETCH)
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

if ~exist('cvx_begin','file')
    error('This test requires <a href="http://cvxr.com/cvx">CVX</a>.');
end
if ~exist(['data/G40.mat'],'file')
    error('Please install G40 dataset from <a href="https://www.cise.ufl.edu/research/sparse/matrices/Gset/G40.html">this link</a> and copy under data folder.');
end

load(['data/G40.mat']);
n = size(Problem.A,1);
L = sparse(diag(Problem.A*ones(n,1)) - Problem.A);
clearvars Problem; 

%IMPORTANT NOTE: I observed that SDPT3 fails with some datasets on
%the maxcut SDP setup below. This is why I changed to the MOSEK solver. 
%If you do not have MOSEK, please check the spectrum of the output of
%CVX. It should not have negative singular values (it is not a problem 
%if they are on the numerical error orders). 
cvx_begin sdp
    cvx_solver mosek
    variable A(n,n) symmetric
    minimize -0.25*trace(A'*L)
    diag(A) == ones(n,1)
    A == semidefinite( n )
cvx_end
clearvars L

[~, singVals, ~] = svd(A);
singVals = diag(singVals);

save(['data/MaxCut']);
