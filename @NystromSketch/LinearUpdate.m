function obj = LinearUpdate( obj, varargin )
%LINEARUPDATE Nystrom Sketch for Fixed-Rank Approximation of a PSD Matrix
%This function implements a linear update routine for the sketches
%described in [TYUC2017Nys].
%   Require: (n x n) dimensinal update matrix (H)
%            scalars (theta1) and (theta2)
%   Ensure:  Modifies sketch (Y) to reflect linear update
%                    A = theta1*A + theta2*H
%            which takes the form
%                    Y = theta1*Y + theta2*H*Omega
%
%   S = S.LinearUpdate(theta1, theta2, H) updates the sketch Y. H can be a
%   matrix of the function handle of a linear operator.
%
%   S = S.LinearUpdate(theta1, theta2, U, V) where U and V are tall
%   matrices updates the sketches Y from the factors H = U*V' in a storage
%   efficient way.
%
%[TYUC2017Nys] J.A. Tropp, A. Yurtsever, M. Udell and V. Cevher. Fixed-Rank
%Approximation of a Positive-Semidefinite Matrix from Streaming Data.
%In Proc. 31st Conference on Neural Information Processing Systems
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

narginchk(2,5);
if nargin == 2
    H = varargin{1};
    theta1 = 0;
    theta2 = 1;
    obj.Y = theta1*obj.Y + theta2*(H*obj.Omega);
elseif nargin == 3
    H = varargin{1};
    theta2 = varargin{2};
    theta1 = 1 - theta2;
    obj.Y = theta1*obj.Y + theta2*(H*obj.Omega);
elseif nargin == 4
    H = varargin{1};
    theta1 = varargin{2};
    theta2 = varargin{3};
    obj.Y = theta1*obj.Y + theta2*(H*obj.Omega);
else
    U       = varargin{1};
    V       = varargin{2};
    theta1  = varargin{3};
    theta2  = varargin{4};
    obj.Y   = theta1*obj.Y + theta2*(U*(V'*obj.Omega));
end

end
