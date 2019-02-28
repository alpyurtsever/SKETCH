function [ Y, S ] = GittensMahoneyApprox( obj, r )
%GITTENSMAHONEYAPPROX Nystrom Sketch for Fixed-Rank Approximation of a PSD Matrix
%This function implements the Fixed-Rank Nystrom Approximation with the standard
%truncation approach (See (2.6) in [TYUC2017Nys]). 
%
%   Require: Rank parameter (1 <= r <= k) 
%   Ensure:  For (q = k + l), returns factors 
%             (n x q) dimensional U with orthonormal columns 
%             (q x q) dimensional non-negative diagonal matrix D 
%            that form a rank-q positive-semidefinite approximation 
%            (Aout = U*D*U') of the sketched matrix. 
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
%Created: April 19, 2017
%Last modified: July 16, 2018 (Modified for SKETCHv1.0)
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
[V, D]  = eig( 0.5*(B+B') , 'vector' );  
[D, I]  = sort(D,'descend');
V       = V(:,I(1:r));
D       = diag(D(1:r));
V       = V(:,1:r);   
D       = D(1:r,1:r);
S       = pinv( subplus(D) );
Y       = obj.Y*V;

if nargout == 1; Y = Y*S*Y'; end

end