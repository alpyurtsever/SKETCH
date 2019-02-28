function [ U, D ] = LowRankSymApprox( obj )
%LOWRANKSYMAPPROX Single-View Low-Rank Symmetric Approximation
%This function implements (5.6) [Algorithm 5] from the main reference.
%   Require: Matrix dimensions to be equal (m = n)
%   Ensure:  For (q = k + l), returns factors 
%             (n x q) dimensional U with orthonormal columns
%             (q x q) dimensional conjugate symmetric matrix D
%            that form a rank-q conjugate symmetric approximation 
%            (Aout = U*D*U') of the sketched matrix.
%
%   Aout = S.LowRankSymApprox() returns the rank-q conjugate symmetric
%   approximation of the square target matrix A for some q <= k. 
%
%   [U, D] = S.SimpleLowRankApprox() returns the factors U and S that
%   form the rank-q conjugate symmetric approximation Aout = U*S*U' for
%   some q <= k. Here, U  is the (m x q) dimensional factor with
%   orthonormal columns, and D is of (q x q) dimensions and conjugate
%   symmetric. 
%
%See our reference paper for the detailed explanation of the sketching
%procedure and the arithmetic, communication and storage costs.
%   
%[TYUC2017] J.A. Tropp, A. Yurtsever, M. Udell and V. Cevher, Practical  
%Sketching Algorithms for Low-Rank Matrix Approximation, SIMAX, 2017.
%
%   NOTICE: THIS FILE IS MODIFIED FROM THE PRACTICALSKETCHING TOOLBOX 
%   by Alp Yurtsever to add SSFT test matrices.
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

if size(obj.Omega,2) ~= size(obj.Upsilon,2)
    error('Matrix should be symmetric (i.e., m = n).')
end
[Q, W] = LowRankApprox(obj);
[U, T] = qr([Q, W'], 0);
k = size(obj.Omega,1);
T1 = T(:,1:k);
T2 = T(:, (k+1):end);
D  = (T1*T2' + T2*T1')/2; 

if nargout == 1; U = U*D*U'; end

end

