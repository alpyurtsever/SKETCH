function [ U, Sigma, V ] = UpaFixedRankApprox( obj, r )
%UPAFIXEDRANKAPPROX Boutisidis' Single-View Fixed-Rank Approximation.
%This function implements Boutisidis' Single-View Fixed-Rank Approximation
%[BWZ2016], using simplifications suggested in [Upa2016]. See our main 
%reference for more details.
%
%   Require: Target rank (r <= k)
%   Ensure:  Returns factors 
%             (m x r) dimensional U with orthonormal columns
%             (r x r) dimensional nonnegative diagonal matrix (Sigma) 
%             (n x r) dimensional V with orthonormal columns
%            that form a rank-r approximation (Aupa = U*Sigma*V') of the 
%            sketched matrix.
%
%   Ar = S.UpaFixedRankApprox(r) returns the rank-r approximation of the
%   target matrix A. 
%
%   [U, Sigma, V] = S.UpaFixedRankApprox(r) returns the factors U, Sigma
%   and V that form the rank-r approximation Aupa = U*Sigma*V'. Here, U
%   and V are (m x r) and (n x r) dimensional factors with orthonormal
%   columns, and Sigma is the (r x r) dimensional matrix. 
%
%See our reference paper [TYUC2019] for the detailed explanation of the
%sketching procedure and the arithmetic and storage costs.
%
%[BWZ2016] C. Boutsidis, D. Woodruff and P. Zhong, Optimal Principal 
%Component Analysis in Distributed and Streaming Models.
%
%[Upa2016] J. Upadhyay, Fast and Space-Optimal Low-Rank Factorization in 
%the Streaming Model with Application in Differential Privacy.
%
%[TYUC2019] J.A. Tropp, A. Yurtsever, M. Udell and V. Cevher.  
%Streaming Low-Rank Matrix Approximation with an Application to
%Scientific Simulation. 
%
%Coded by: Alp Yurtsever
%Ecole Polytechnique Federale de Lausanne, Switzerland.
%Laboratory for Information and Inference Systems, LIONS.
%contact: alp.yurtsever@epfl.ch
%Created: December 30, 2016 
%Last modified: July 16, 2018 (Modified from PRACTICALSKETCHINGv1.0)
%
% SKETCHv1.1
% Copyright (C) 2018 Laboratory for Information and Inference Systems
% (LIONS), Ecole Polytechnique Federale de Lausanne, Switzerland.
% This code is a part of SKETCH toolbox. 
% Please read COPYRIGHT before using this file.

if r > min(size(obj.Omega, 2), size(obj.Psi,2))
    error('Target rank ''r'' cannot be greater than matrix dimensions.')
end

% Computing the approximation
[P, ~]   = qr(obj.X',0);
[Q, ~]   = qr(obj.Y,0); 
[U1, T1]  = qr(obj.Phi*Q,0);
[U2, T2]  = qr(obj.Psi*P,0);

[Uupa, Sigma, Vupa] = svd(U1'*obj.Z*U2);
Uupa = Uupa(:,1:r);
Sigma = Sigma(1:r,1:r);
Vupa = Vupa(:,1:r);

W = (T1\Uupa)*Sigma*(Vupa'/T2');
[Uw,Sigma,Vw] = svd(W);

U = Q*Uw;
V = P*Vw;

if nargout == 1; U = U*Sigma*V'; end

end

