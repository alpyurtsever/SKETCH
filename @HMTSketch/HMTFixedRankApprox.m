function [ U, Sigma, V ] = HMTFixedRankApprox( obj, r )
%HMTFIXEDRANKAPPROX Single-View Fixed-Rank Approximation
%This function implements the fixed rank approximation from [HMT2011].
%
%   Require: Target rank (r <= k)
%   Ensure:  Returns factors 
%             (m x r) dimensional U with orthonormal columns
%             (n x r) dimensional V with orthonormal columns
%             (r x r) dimensional nonnegative diagonal matrix (Sigma) 
%            that form a rank-r approximation (Aout = U*Sigma*V') of the 
%            sketched matrix.
%
%   Ar = S.HMTFixedRankApprox(r) returns the rank-r approximation of the
%   target matrix A. 
%
%   [U, Sigma, V] = S.HMTFixedRankApprox(r) returns the factors U, Sigma
%   and V that form the rank-r apptroximation Ar = U*Sigma*V'. Here, U
%   and V are (m x r) and (n x r) dimensional factors with orthonormal
%   columns, and Sigma is the (r x r) dimensional nonnegative diagonal
%   matrix. 
%   
%[HMT2011] N. Halko, P.G. Martinsson and J.A. Tropp. Finding Structure
%with Randomness: Probabilistic algorithms for constructing approximate 
%matrix decompositions,
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
%SKETCHv1.1
%Copyright (C) 2018 Laboratory for Information and Inference Systems
%(LIONS), Ecole Polytechnique Federale de Lausanne, Switzerland.
%This code is a part of SKETCH toolbox. 
%Please read COPYRIGHT before using this file.

if r > min(size(obj.Omega, 1), size(obj.OmegaTilde,1))
    error('Target rank ''r'' cannot be greater than matrix dimensions.')
end

[Q,~,~] = svd(obj.Y,0);             % Implement HMT11, Remark 5.4
Q = Q(:, 1:r);

[Qtilde,~,~] = svd(obj.Ytilde,0);   % Implement HMT11, Remark 5.4
Qtilde = Qtilde(:, 1:r);

B1 = (Q'*obj.Y) / (Qtilde'*obj.Omega);
B2 = ((Q'*obj.OmegaTilde)') \ ((Qtilde'*obj.Ytilde)');

X = 0.5*(B1+B2);                    % LS solution to HMT (5.14), (5.15)

[Ux,Sigma,Vx] = svd(X);
U = Q*Ux;
V = Qtilde*Vx;

if nargout == 1; U = U*Sigma*V'; end

end