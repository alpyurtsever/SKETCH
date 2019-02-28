function [ U, D ] = FixedRankPSDApprox( obj, r )
%FIXEDRANKPSDAPPROX Single-View Fixed-Rank PSD Approximation
%This function implements (6.6) [Algorithm 9] from the main reference.
%   Require: Matrix dimensions to be equal (m = n)
%            Target rank (r <= k)
%   Ensure:  Returns factors 
%             (n x r) dimensional U with orthonormal columns
%             (r x r) dimensional nonnegative diagonal matrix D
%            that form a rank-r positive-semidefinite (PSD) approximation 
%            (Aout = U*D*U') of the sketched matrix.
%
%   Aout = S.FixedRankPSDApprox(r) returns the rank-r positive
%   semi-definite approximation of the square target matrix A.
%
%   [U, D] = S.FixedRankPSDApprox(r) returns the factors U and D that
%   form the rank-r positive semi-definite approximation Aout = U*D*U'.
%   Here, U  is the (m x q) dimensional factor with orthonormal
%   columns, and D is of (q x q) dimensions and non-negative diagonal
%   matrix.
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

if r > min(size(obj.Omega,2), size(obj.Upsilon,2))
    error('Target rank ''r'' cannot be greater than matrix dimensions.')
end

[U, S] = LowRankSymApprox( obj );

% [V, D] = eigs(S, r, 'LA'); % This does not work when S is complex valued
% even if it is conugate symmetric. We can always use 'LR', which works
% fine but gives a warning when S is real symmetric. We can check these
% conditions and use 'la' and 'lr' accordingly. But the code below is
% numerically more stable! Note that S is a very small matrix!
[V, D] = eig(S);
D = diag(D);
D = real(D); % D comes complex due to numerical errors!
[D,I] = sort(D,1,'descend');
V = V(:,I(1:r));
D = diag(D(1:r));
U = U*V;
D = subplus(D);

if nargout == 1; U = U*D*U'; end

end


