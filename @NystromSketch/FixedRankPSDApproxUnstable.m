function [ U, Delta ] = FixedRankPSDApproxUnstable( obj, r )
%FIXEDRANKPSDAPPROX Fixed-Rank Approximation from Nystrom Sketch.
%This function implements Fixed-Rank PSD Approximation approach proposed in  
%[LLSSKT2017, Eqn.(13)]. 
%
%   Require: Rank parameter (1 <= r <= k) 
%   Ensure:  For (q = k + l), returns factors 
%             (n x q) dimensional U with orthonormal columns 
%             (q x q) dimensional non-negative diagonal matrix Delta 
%            that form a rank-q positive-semidefinite approximation 
%            (Aout = U*Delta*U') of the sketched matrix. 
%
%[LLSSKT2017] H.Li, G.C.Linderman, A.Szlam, K.P.Stanton, Y.Kluger, M. Tygert
%Algorithm 971: An implementation of a randomized algorithm for principal
%component analysis. ACM Trans. Math. Softw., 43(3):28:1???28:14, Jan. 2017
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

if (r < 1) || (r > obj.k)
    error('Target rank ''r'' must be a positive integer less than the sketch size parameter ''k''.')
end

B = obj.Omega'*obj.Y;
[V, Delta]  = eig(0.5*(B+B'));
S           =  V*pinv(sqrt(subplus(Delta)))*V';
[U, Sigma, ~] = svd(obj.Y*S,'econ'); U = U(:,1:r); Sigma = Sigma(1:r,1:r);
Delta       = Sigma.^2;

if nargout == 1; U = U*Delta*U'; end

end
