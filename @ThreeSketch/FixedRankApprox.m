function [ U, Sigma, V ] = FixedRankApprox( obj, r )
%FIXEDRANKPSDAPPROX Single-View Fixed-Rank Approximation
%This function implements the fixed-rank approximation method explained in 
%[TYUC2019]. See our main reference for more details.
%
%   Require: Target rank (r <= k)
%   Ensure:  Returns factors 
%             (m x r) dimensional U with orthonormal columns
%             (r x r) dimensional nonnegative diagonal matrix D
%             (n x r) dimensional V with orthonormal columns
%            that form a rank-r positive-semidefinite (PSD) approximation 
%            (Aout = U*D*U') of the sketched matrix.
%
%   Aout = S.FixedRankPSDApprox(r) returns the rank-r approximation
%   of the target matrix A.
%
%   [U, D, V] = S.FixedRankPSDApprox(r) returns the factors U, D and V 
%   that form the rank-r approximation Aout = U*D*V'.
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

[Q, W, P] = LowRankApprox( obj );
try
    [U,Sigma,V] = svd(W);
    U = U(:,1:r);
    V = V(:,1:r);
    Sigma = Sigma(1:r,1:r);
catch
    [U,Sigma,V] = svds(W,r);
end
U = Q*U;
V = P*V;

if nargout == 1; U = U*Sigma*V'; end

end

