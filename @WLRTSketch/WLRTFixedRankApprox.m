function [ U, Sigma, V ] = WLRTFixedRankApprox( obj, r )
%WLRTFIXEDRANKAPPROX Single-View Fixed-Rank Approximation
%This function implements the algorithm in Section 5.2 from [WLRT2008].
%
%   Require: Target rank (r <= k)
%   Ensure:  Returns factors 
%             (m x r) dimensional Q with orthonormal columns
%             (n x r) dimensional V with orthonormal columns
%             (r x r) dimensional nonnegative diagonal matrix (Sigma) 
%            that form a rank-r approximation (Aout = U*Sigma*V') of the 
%            sketched matrix.
%
%   Ar = S.WLRTFixedRankApprox(r) returns the rank-r approximation of the
%   target matrix A. 
%
%   [Q, Sigma, V] = S.WLRTFixedRankApprox(r) returns the factors Q, Sigma
%   and V that form the rank-r apptroximation Ar = U*Sigma*V'. Here, U
%   and V are (m x r) and (n x r) dimensional factors with orthonormal
%   columns, and Sigma is the (r x r) dimensional nonnegative diagonal
%   matrix. 
%
%[WLRT2008] F. Woolfe, E. Liberty, V. Rokhlin, M. Tygert. A Fast
%Randomized Algorithm for the Approximation of Matrices.
%
%[TYUC2019] J.A. Tropp, A. Yurtsever, M. Udell and V. Cevher.  
%Streaming Low-Rank Matrix Approximation with an Application to
%Scientific Simulation. 
%
%Coded by: Alp Yurtsever
%Ecole Polytechnique Federale de Lausanne, Switzerland.
%Laboratory for Information and Inference Systems, LIONS.
%contact: alp.yurtsever@epfl.ch
%Created: January 18, 2019
%Last modified: January 22, 2019 by JAT
%
% SKETCHv1.1
% Copyright (C) 2018 Laboratory for Information and Inference Systems
% (LIONS), Ecole Polytechnique Federale de Lausanne, Switzerland.
% This code is a part of SKETCH toolbox. 
% Please read COPYRIGHT before using this file.

if r > min(size(obj.R, 2), size(obj.Rtilde,2))
    error('Target rank ''r'' cannot be greater than matrix dimensions.')
end

[Q,~,~] = svd(obj.Y',0);
Q = Q(:,1:r);
[P,~,~] = svd(obj.Ytilde',0);
P = P(:,1:r);
W = obj.R*P;
B = obj.Y*Q;
% X = W\B;

[WQ,WR] = qr(W,0);      % JAT: Changed notation for clarity
X = WR\(WQ'*B);         

[Ux,Sigma,Vx] = svd(X);
U = P*Ux;
V = Q*Vx;

if nargout == 1; U = U*Sigma*V'; end

end

