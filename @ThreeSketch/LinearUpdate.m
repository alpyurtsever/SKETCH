function LinearUpdate( obj, varargin )
%LINEARUPDATE Single-View Sketch: Linear Update
%This function implements the linear updates explained in [TYUC2019]. 
%See our main reference for more details.
%
%   Require: (m x n) dimensinal update matrix (H)
%            scalars (theta) and (tau)
%   Ensure:  Modifies sketch (Y,W) to reflect linear update
%                    A = theta*A + eta*H
%            which takes the form
% 					 X = theta*X + tau*Upsilon*H
%                    Y = theta*Y + tau*H*Omega'
%                    Z = theta*Z + tau*Psi*H*phi'
%
%   S = LinearUpdate(H) updates the sketches X, Y and Z by choosing
%	(tau = 1) and (theta = 0).
%
%   S = LinearUpdate(H, tau) updates the sketches X, Y and Z by
%	choosing (theta = 1 - tau).
%
%   S = LinearUpdate(H, theta, tau) updates the sketches X, Y and Z.
%
%   S = LinearUpdate(U, V, theta, tau) where U and V are tall
%   matrices updates the sketches X, Y and Z from the factors H = U*V'
%   in an efficient way.
%
% [TYUC2019] J.A. Tropp, A. Yurtsever, M. Udell and V. Cevher.  
% Streaming Low-Rank Matrix Approximation with an Application to
% Scientific Simulation. 
%
%Coded by: Alp Yurtsever
%Ecole Polytechnique Federale de Lausanne, Switzerland.
%Laboratory for Information and Inference Systems, LIONS.
%contact: alp.yurtsever@epfl.ch
%Created: March 01, 2018
%Last modified: February 11, 2019
%
% SKETCHv1.1
% Copyright (C) 2018 Laboratory for Information and Inference Systems
% (LIONS), Ecole Polytechnique Federale de Lausanne, Switzerland.
% This code is a part of SKETCH toolbox. 
% Please read COPYRIGHT before using this file.

narginchk(2,5);

m = size(obj.Y,1); %#ok
n = size(obj.X,2);

if nargin == 2
    H = varargin{1};
    eta = 0;
    nu = 1;
    obj.X = eta*obj.X + nu*(obj.Upsilon*H);
    obj.Y = eta*obj.Y + nu*(H*obj.Omega');
    obj.Z = eta*obj.Z + nu*(obj.Phi*H*obj.Psi');
    if ~isempty(obj.Theta)
        obj.W = eta*obj.W + nu*(obj.Theta*H);
    end
elseif nargin == 3
    H     = varargin{1};
    nu    = varargin{2};
    eta   = 1 - nu;
    obj.X = eta*obj.X + nu*(obj.Upsilon*H);
    obj.Y = eta*obj.Y + nu*(H*obj.Omega');
    obj.Z = eta*obj.Z + nu*(obj.Phi*H*obj.Psi');
    if ~isempty(obj.Theta)
        obj.W = eta*obj.W + nu*(obj.Theta*H);
    end
elseif nargin == 4
    H      = varargin{1};
    eta    = varargin{2};
    nu     = varargin{3};
    obj.X = eta*obj.X + nu*(obj.Upsilon*H);
    obj.Y = eta*obj.Y + nu*(H*obj.Omega');
    obj.Z = eta*obj.Z + nu*(obj.Phi*H*obj.Psi');
    if ~isempty(obj.Theta)
        obj.W = eta*obj.W + nu*(obj.Theta*H);
    end
else
    Hforw  = varargin{1};
    Hback  = varargin{2};
    eta    = varargin{3};
    nu     = varargin{4};
    obj.X = eta*obj.X + nu*((obj.Upsilon*Hforw)*Hback');
    obj.Y = eta*obj.Y + nu*(Hforw*(Hback'*obj.Omega'));
    obj.Z = eta*obj.Z + nu*((obj.Phi*Hforw)*(Hback'*obj.Psi'));
    if ~isempty(obj.Theta)
        obj.W = eta*obj.W + nu*((obj.Theta*Hforw)*Hback');
    end
end


end

