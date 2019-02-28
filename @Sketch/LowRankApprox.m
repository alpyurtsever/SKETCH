function [ Q, W ] = LowRankApprox( obj )
%LOWRANKAPPROX Single-View Low-Rank Approximation
%This function implements (4.3) [Algorithm 4] from the main reference.
%   Ensure:  Returns factors 
%             (m x k) dimensional Q with orthonormal columns
%             (k x n) dimensional W
%            that form a rank-k approximation (Aout = Q*W) of the sketched
%            matrix.
%
%   Aq = S.LowRankApprox() returns a rank-q approximation of the target
%   matrix A for some q <= k. 
%
%   [Q, W] = S.LowRankApprox() returns the factors Q and W that
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

[Q, ~] = qr(obj.Y,0);
[U,T] = qr(obj.Upsilon*Q,0);
W = T\(U'*obj.X);

if nargout == 1; Q = Q*W; end

end

