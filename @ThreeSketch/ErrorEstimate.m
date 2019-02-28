function [ err2 ] = ErrorEstimate( obj, varargin )
%ERRORSKETCH implements error estimation procedure for sketch with three
%components. This function returns the error estimate of the input matrix 
%Ahat, as described in Algorithm 6.1 in [TYUC2019]. See our main reference 
%for more details.
%
%   Require: Approximation Ahat (or in the factored form as Q*W*P)
%   Ensure:  Returns the error estimate. 
%
%   err = S.ErrorEstimate(Ahat) returns the error estimate for Ahat.
%
%   err = S.ErrorEstimate(Q,W,P) returns the error estimate for Q*W*P'.
%
%   err0 = S.ErrorEstimate(0) returns the estimate of norm(A,'fro')^2.
%   
% [TYUC2019] J.A. Tropp, A. Yurtsever, M. Udell and V. Cevher.  
% Streaming Low-Rank Matrix Approximation with an Application to
% Scientific Simulation. 
%
%Coded by: Alp Yurtsever
%Ecole Polytechnique Federale de Lausanne, Switzerland.
%Laboratory for Information and Inference Systems, LIONS.
%contact: alp.yurtsever@epfl.ch
%Created: February 05, 2019
%Last modified: February 22, 2019
%
% SKETCHv1.1
% Copyright (C) 2018 Laboratory for Information and Inference Systems
% (LIONS), Ecole Polytechnique Federale de Lausanne, Switzerland.
% This code is a part of SKETCH toolbox. 
% Please read COPYRIGHT before using this file.

if isreal(obj.Theta)
    beta = 1;
else
    beta = 2;
end

q = size(obj.Theta,1);

if nargin == 2
    Aout = varargin{1};
    if Aout == 0
        err2 = (1/sqrt(beta*q)) * norm(obj.W,'fro');
    else
        err2 = (1/sqrt(beta*q)) * norm(obj.W-obj.Theta*Aout,'fro');
    end
elseif nargin == 4
    Q = varargin{1};
    W = varargin{2};
    P = varargin{3};
    err2 = (1/sqrt(beta*q)) * norm(obj.W-((obj.Theta*Q)*W)*P','fro');
end

end

