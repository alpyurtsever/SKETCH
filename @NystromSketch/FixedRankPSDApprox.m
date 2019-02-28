function [ U, Delta ] = FixedRankPSDApprox( obj, r )
%FIXEDRANKPSDAPPROX Fixed-Rank Approximation from Nystrom Sketch.
%This function implements Algorithm 3 from the main reference [TYUC2017Nys]. 
%
%   Require: Rank parameter (1 <= r <= k) 
%   Ensure:  For (q = k + l), returns factors 
%             (n x q) dimensional U with orthonormal columns 
%             (q x q) dimensional non-negative diagonal matrix Delta 
%            that form a rank-q positive-semidefinite approximation 
%            (Aout = U*Delta*U') of the sketched matrix. 
%   
%[TYUC2017Nys] J.A. Tropp, A. Yurtsever, M. Udell and V. Cevher. Fixed-
%Rank Approximation of a Positive-Semidefinite Matrix from Streaming 
%Data. In Proc. 31st Conference on Neural Information Processing Systems
%(NeurIPS), Long Beach, CA, USA, December 2017.
%
%Coded by: Alp Yurtsever
%Ecole Polytechnique Federale de Lausanne, Switzerland.
%Laboratory for Information and Inference Systems, LIONS.
%contact: alp.yurtsever@epfl.ch
%Created: April 12, 2017
%Last modified: July 16, 2018 (Modified for SKETCHv1.0, by Alp Yurtsever)
%
%NysSKETCHv1.0
%Copyright (C) 2017 Laboratory for Information and Inference Systems
%(LIONS), Ecole Polytechnique Federale de Lausanne, Switzerland.
%This code is a part of Nys???SKETCH toolbox. 
%Please read COPYRIGHT before using this file.
%This code is a part of SKETCH toolbox. 
%Please read COPYRIGHT before using this file.

if (r < 1) || (r > size(obj.Omega,2))
    error('Target rank ''r'' must be a positive integer less than the sketch size parameter ''k''.')
end

Y = obj.Y;
nu = eps*norm(Y);

Y = Y + nu*obj.Omega;
B = obj.Omega' * Y;
    
B = 0.5*(B+B');
C = chol(B);
[U, Sigma, ~] = svd( Y / C, 'econ' );
U = U(:, 1:r); 
Sigma = Sigma( 1:r, 1:r );
Delta = subplus( Sigma^2 - nu*eye(r) );

if nargout == 1; U = U*Delta*U'; end

end

