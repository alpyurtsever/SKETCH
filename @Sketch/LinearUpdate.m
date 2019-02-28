function obj = LinearUpdate( obj, varargin )
%LINEARUPDATE Single-View Sketch: Linear Update
%This function implements (3.6) from the main reference.
%   Require: (m x n) dimensinal update matrix (H)
%            scalars (theta) and (eta)
%   Ensure:  Modifies sketch (Y,W) to reflect linear update
%                    A = theta*A + eta*H
%            which takes the form
%                    Y = theta*Y + eta*H*Omega'
%                    W = theta*W + eta*Upsilon*H
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
%   updates the sketches Y and W as Y = theta*Y + eta*H*Omega' and 
%   W = theta*W + eta*Upsilon*H.
%
%See our reference paper for the detailed explanation of the sketching
%procedure and the arithmetic, communication and storage costs.
%   
%[TYUC2017] J.A. Tropp, A. Yurtsever, M. Udell and V. Cevher, Practical  
%Sketching Algorithms for Low-Rank Matrix Approximation, 2017.
%
%Coded by: Alp Yurtsever
%Ecole Polytechnique Federale de Lausanne, Switzerland.
%Laboratory for Information and Inference Systems, LIONS.
%contact: alp.yurtsever@epfl.ch
%Created: June 28, 2016
%Last modified: July 16, 2018 (Modified for SKETCHv1.0, by Alp Yurtsever)
%
%This code is a part of PRACTICALSKETCHING toolbox. 
%Please read COPYRIGHT before using this file.
%This code is a part of SKETCH toolbox. 
%Please read COPYRIGHT before using this file.

narginchk(2,5);

if nargin == 2
    H = varargin{1};
    theta = 0;
    eta = 1;
    obj.Y = theta*obj.Y + eta*(H*obj.Omega');
    obj.X = theta*obj.X + eta*(obj.Upsilon*H);
elseif nargin == 3
    H = varargin{1};
    eta = varargin{2};
    theta = 1 - eta;
    obj.Y = theta*obj.Y + eta*(H*obj.Omega');
    obj.X = theta*obj.X + eta*(obj.Upsilon*H);
elseif nargin == 4
    H = varargin{1};
    theta = varargin{2};
    eta = varargin{3};
    obj.Y = theta*obj.Y + eta*(H*obj.Omega');
    obj.X = theta*obj.X + eta*(obj.Upsilon*H);
else
    Hforw   = varargin{1};
    Hback   = varargin{2};
    theta   = varargin{3};
    eta     = varargin{4};
    obj.Y = theta*obj.Y + eta*(Hforw*(Hback'*obj.Omega'));
    obj.X = theta*obj.X + eta*((obj.Upsilon*Hforw)*Hback');
end

end

