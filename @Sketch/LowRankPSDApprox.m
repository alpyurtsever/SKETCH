function [ U, D ] = LowRankPSDApprox( obj )
%LOWRANKPSDAPPROX Single-View Low-Rank Positive-Semidefinite Approximation
%This function implements (5.7) [Algoritm 6] from the main reference.
%   Require: Matrix dimensions to be equal (m = n)
%   Ensure:  For (q = k + l), returns factors 
%             (n x q) dimensional U with orthonormal columns
%             (q x q) dimensional non-negative diagonal matrix D
%            that form a rank-q positive-semidefinite approximation 
%            (Aout = U*D*U') of the sketched matrix.
%
%   Aout = S.LowRankPSDApprox() returns the rank-q positive
%   semi-definite approximation of the square target matrix A for some
%   q <= k. 
%
%   [U, D] = S.LowRankPSDApprox() returns the factors U and D that
%   form the rank-q positive semi-definite approximation Aout = U*D*U'
%   for some q <= k. Here, U  is the (m x q) dimensional factor with
%   orthonormal columns, and D is of (q x q) dimensions and
%   non-negative diagonal matrix. 
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

if size(obj.Omega, 2) ~= size(obj.Upsilon,2)
    error('Matrix should be symmetric (i.e., m = n).')
end

[U, S] = LowRankSymApprox(obj);
[V, D] = eig(S);
%D = real(D); % D comes complex due to numerical errors!
U = U*V;
D = subplus(D);

if nargout == 1; U = U*D*U'; end

end

