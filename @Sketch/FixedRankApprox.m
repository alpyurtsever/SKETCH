function [ Q, Sigma, V ] = FixedRankApprox( obj, r )
%FIXEDRANKAPPROX Single-View Fixed-Rank Approximation
%This function implements (6.3) [Algorithm 7] from the main reference.
%   Require: Target rank (r <= k)
%   Ensure:  Returns factors 
%             (m x r) dimensional Q with orthonormal columns
%             (n x r) dimensional V with orthonormal columns
%             (r x r) dimensional nonnegative diagonal matrix (Sigma) 
%            that form a rank-r approximation (Aout = Q*Sigma*V') of the 
%            sketched matrix.
%
%   Ar = S.FixedRankApprox(r) returns the rank-r approximation of the
%   target matrix A. 
%
%   [Q, Sigma, V] = S.FixedRankApprox(r) returns the factors Q, Sigma
%   and V that form the rank-r apptroximation Ar = U*Sigma*V'. Here, U
%   and V are (m x r) and (n x r) dimensional factors with orthonormal
%   columns, and Sigma is the (r x r) dimensional nonnegative diagonal
%   matrix. 
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

if r > min(size(obj.Omega, 2), size(obj.Upsilon,2))
    error('Target rank ''r'' cannot be greater than matrix dimensions.')
end

[Q, W] = LowRankApprox(obj);
if min(size(W)) > 5000
    [U, Sigma, V] = svds(W, r);
else
    [U, Sigma, V] = svd(W, 'econ');
    U = U(:,1:r);
    Sigma = Sigma(1:r,1:r);
    V = V(:,1:r);
end
Q = Q*U;

if nargout == 1; Q = Q*Sigma*V'; end

end

