function [ U, D ] = FixedRankSymApprox( obj, r )
%FIXEDRANKSYMAPPROX Single-View Fixed-Rank Symmetric Approximation
%This function implements (6.5) [Algorithm 8] from the main reference.
%   Require: Matrix dimensions to be equal (m = n)
%            Target rank (r <= k)
%   Ensure:  Returns factors 
%             (n x r) dimensional U with orthonormal columns
%             (r x r) dimensional diagonal conjugate symmetric matrix D
%            that form a rank-r conjugate symmetric approximation 
%            (Aout = U*D*U') of the sketched matrix.
%
%   Aout = S.FixedRankSymApprox(r) returns the rank-r conjugate
%   symmetric approximation of the square target matrix A. 
%
%   [U, D] = S.FixedRankSymApprox(r) returns the factors U and S that
%   form the rank-r conjugate symmetric approximation Aout = U*D*U'.
%   Here, U  is the (m x q) dimensional factor with orthonormal
%   columns, and D is of (q x q) dimensions and conjugate symmetric.
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
    error('Target matrix should be square (i.e., m = n).')
end
if r > min(size(obj.Omega, 2), size(obj.Upsilon,2))
    error('Target rank ''r'' cannot be greater than matrix dimensions.')
end

[U, S] = LowRankSymApprox( obj );

% [V, D] = eigs(S, r, 'LM'); % This line is correct, but the code below is
% numerically more stable! Note that D is a very small matrix!
[V, D] = eig(S);
D = diag(D);
D = real(D); % D may come complex due to numerical errors!
[D, I] = sort(abs(D),1,'descend');
V = V(:,I(1:r));
D = diag(D(1:r));
U = U*V;

if nargout == 1; U = U*D*U'; end

end

