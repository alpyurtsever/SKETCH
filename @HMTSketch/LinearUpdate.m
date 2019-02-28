function obj = LinearUpdate( obj, varargin )
%LINEARUPDATE HMTSketch: Linear Update
%This function implements linear updates for HMT-type sketch.
%   Require: (m x n) dimensinal update matrix (H)
%            scalars (theta) and (eta)
%   Ensure:  Modifies sketch (Y,W) to reflect linear update
%                    A = theta*A + eta*H
%            which takes the form
%                    Y = theta*Y + eta*H*Omega
%                    Ytilde = theta*W + eta*H'*OmegaTilde
%
%   S = S.LinearUpdate(H, theta, eta) updates the sketches Y and W as 
%   Y = theta*Y + eta*H*Omega and Ytilde = theta*W + eta*H'*OmegaTilde.
%
%   S = S.LinearUpdate(U, V, theta, eta) where Hf and Hb are tall
%   matrices updates the sketches Y and Ytilde from the factors H = U*V' 
%   in a storage efficient way.
%   
%   S = S.LinearUpdate(Hforw, Hback, theta, eta) where Hforw and Hback 
%   are the function handles to compute H*x = Hforw(x) and x*H = Hback(x),
%   updates the sketches Y and Ytilde in an efficient way.
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
%Last modified: January 22, 2019
%
% SKETCHv1.1
% Copyright (C) 2018 Laboratory for Information and Inference Systems
% (LIONS), Ecole Polytechnique Federale de Lausanne, Switzerland.
% This code is a part of SKETCH toolbox. 
% Please read COPYRIGHT before using this file.

narginchk(2,5);

if nargin == 2
    H = varargin{1};
    theta = 0;
    eta = 1;
    obj.Y = theta*obj.Y + eta*(H*obj.Omega);
    obj.Ytilde = theta*obj.Ytilde + eta*(H'*obj.OmegaTilde);
elseif nargin == 3
    H = varargin{1};
    eta = varargin{2};
    theta = 1 - eta;
    obj.Y = theta*obj.Y + eta*(H*obj.Omega);
    obj.Ytilde = theta*obj.Ytilde + eta*(H'*obj.OmegaTilde);
elseif nargin == 4
    H = varargin{1};
    theta = varargin{2};
    eta = varargin{3};
    obj.Y = theta*obj.Y + eta*(H*obj.Omega);
    obj.Ytilde = theta*obj.Ytilde + eta*(H'*obj.OmegaTilde);
else
    Hforw   = varargin{1};
    Hback   = varargin{2};
    theta   = varargin{3};
    eta     = varargin{4};
    obj.Y = theta*obj.Y + eta*(Hforw*(Hback'*obj.Omega));
    obj.Ytilde = theta*obj.Ytilde + eta*(Hback*(Hforw'*obj.OmegaTilde));
end

end

