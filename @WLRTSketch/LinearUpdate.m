function obj = LinearUpdate( obj, varargin )
%LINEARUPDATE WLRTSketch: Linear Update
%This function implements linear updates for WLRT-type sketch.
%   Require: (m x n) dimensinal update matrix (H)
%            scalars (theta) and (eta)
%   Ensure:  Modifies sketch (Y,W) to reflect linear update
%                    A = theta*A + eta*H
%            which takes the form
%                    Y = theta*Y + eta*R*H
%                    W = theta*W + eta*Rtilde*H'
%
%   S = S.LinearUpdate(H, theta, eta) updates the sketches Y and W as 
%   Y = theta*Y + eta*H*Omega' and W = theta*W + eta*Upsilon*H.
%
%   S = S.LinearUpdate(U, V, theta, eta) where Hf and Hb are tall
%   matrices updates the sketches Y and W from the factors H = U*V' 
%   in a storage efficient way.
%   
%   S = S.LinearUpdate(Hforw, Hback, theta, eta) where Hforw and Hback 
%   are the function handles to compute H*x = Hforw(x) and x*H = Hback(x),
%   updates the sketches Y and Ytilde in an efficient way.
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
%Created: June 28, 2016
%Last modified: February 21, 2019
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
    obj.Y = theta*obj.Y + eta*(obj.R*H);
    obj.Ytilde = theta*obj.Ytilde + eta*(obj.Rtilde*H');
elseif nargin == 3
    H = varargin{1};
    eta = varargin{2};
    theta = 1 - eta;
    obj.Y = theta*obj.Y + eta*(obj.R*H);
    obj.Ytilde = theta*obj.Ytilde + eta*(obj.Rtilde*H');
elseif nargin == 4
    H = varargin{1};
    theta = varargin{2};
    eta = varargin{3};
    obj.Y = theta*obj.Y + eta*(obj.R*H);
    obj.Ytilde = theta*obj.Ytilde + eta*(obj.Rtilde*H');
else
    Hforw   = varargin{1};
    Hback   = varargin{2};
    theta   = varargin{3};
    eta     = varargin{4};
    obj.Y = theta*obj.Y + eta*((obj.R*Hforw)*Hback');
    obj.Ytilde = theta*obj.Ytilde + eta*((obj.Rtilde*Hback)*Hforw');
end

end

