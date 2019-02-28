function [ Q, W ] = SimpleLowRankApprox( obj )
%SIMPLELOWRANKAPPROX Simplest Single-View Low-Rank Approximation
%This function implements Algorithm 3 from [TYUC2017].
%   Ensure:  For some (q <= k), returns factors 
%             (m x q) dimensional Q with orthonormal columns
%             (q x n) dimensional W
%            that form a rank-q approximation (Aout = Q*W) of the sketched
%            matrix.
%
%   Aq = S.SimpleLowRankApprox() returns a rank-q approximation of
%   the target matrix A for some q <= k.
%
%   [Q, W] = S.SimpleLowRankApprox() returns the factors Q and W that
%   form the rank-q approximation Aq = Q*W for some q <= k. Here, Q is
%   the (m x q) dimensional factor with orthonormal columns, and W is
%   of (q x n) dimensions.
%
%See our reference paper for the detailed explanation of the sketching
%procedure and the arithmetic, communication and storage costs.
%   
%[TYUC2017] J.A. Tropp, A. Yurtsever, M. Udell and V. Cevher, Practical  
%Sketching Algorithms for Low-Rank Matrix Approximation, SIMAX, 2017.
%
%Coded by: Alp Yurtsever
%Ecole Polytechnique Federale de Lausanne, Switzerland.
%Laboratory for Information and Inference Systems, LIONS.
%contact: alp.yurtsever@epfl.ch
%Created: June 28, 2016
%Last modified: July 16, 2018 (Modified for SKETCHv1.0, by Alp Yurtsever)
%
%PRACTICALSKETCHING-v1.0
%Copyright (C) 2017 Laboratory for Information and Inference Systems
%(LIONS), Ecole Polytechnique Federale de Lausanne, Switzerland.
%
%This code is a part of PRACTICALSKETCHING toolbox. 
%Please read COPYRIGHT before using this file.
%This code is a part of SKETCH toolbox. 
%Please read COPYRIGHT before using this file.

Q = orth(obj.Y);
W = (obj.Upsilon*Q)\obj.X;

if nargout == 1; Q = Q*W; end

end

