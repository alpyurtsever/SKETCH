function [ Q, W, P ] = LowRankApprox( obj )
%LOWRANKAPPROX Single-View Low-Rank Approximation
%This function implements the low-rank approximation method from the 
%sketch with three components, explained in TYUC2019]. See our main 
%reference for more details.
%
%   Ensure:  Returns factors 
%             (m x k) dimensional Q with orthonormal columns
%             (k x k) dimensional W
%             (n x k) dimensional P with orthonormal columns
%			 that form a low-rank approximation (Aout = Q*W*P').
%
%   Aq = S.LowRankApprox() returns a rank-q approximation of the target
%   matrix A for some q <= k. 
%
%   [Q, W, P] = S.LowRankApprox() returns the factors Q, W and P that
%   form the rank-q approximation (Aout = Q*W*P') some q <= k.
%
% [TYUC2019] J.A. Tropp, A. Yurtsever, M. Udell and V. Cevher.  
% Streaming Low-Rank Matrix Approximation with an Application to
% Scientific Simulation. 
%
%Coded by: Alp Yurtsever
%Ecole Polytechnique Federale de Lausanne, Switzerland.
%Laboratory for Information and Inference Systems, LIONS.
%contact: alp.yurtsever@epfl.ch
%Created: March 01, 2018
%Last modified: July 16, 2018
%
% SKETCHv1.1
% Copyright (C) 2018 Laboratory for Information and Inference Systems
% (LIONS), Ecole Polytechnique Federale de Lausanne, Switzerland.
% This code is a part of SKETCH toolbox. 
% Please read COPYRIGHT before using this file.

[Q, ~] = qr(obj.Y,0);
[P, ~] = qr(obj.X',0);
[U1,T1] = qr(obj.Phi*Q,0);
[U2,T2] = qr(obj.Psi*P,0);
W = T1\(U1'*obj.Z*U2)/T2';

if nargout == 1; Q = Q*W*P'; end

end

