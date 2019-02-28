function [ U, S, V ] = updateThinSVD( U1, S1, V1, A, B )
% FUNCTION:   [U, S, V] = updateThinSVD(...)
% PURPOSE:    Modification of given SVD under low-rank perturbations.
%             [U, S, V] = svd(U1*S1*V1' + A*B');
% REFERENCE:  [Bra2006] Brand,M., Fast low-rank modification of thin 
%             singular value decomposition, Linear Algebra and its 
%             Applications 415 (2006) 20-30, ELSEVIER, 
%             doi:10.1016/j.laa.2005.07.021
%
%   This function implements the thin SVD update method from [Bra2006].
%   Implemented to be used in the large scale experiments of [TYUC2017]
%
%   [TYUC2017] J.A. Tropp, A. Yurtsever, M. Udell and V. Cevher, Practical  
%   Sketching Algorithms for Low-Rank Matrix Approximation, SIMAX, 2017.
%
%   Coded by: Alp Yurtsever
%   Ecole Polytechnique Federale de Lausanne, Switzerland.
%   Laboratory for Information and Inference Systems, LIONS.
%   contact: alp.yurtsever@epfl.ch
%   Created: July 29, 2016
%   Last modified: August 29, 2016
%
% PRACTICALSKETCHINGv1.0
% Copyright (C) 2017 Laboratory for Information and Inference Systems (LIONS)
%		    Ecole Polytechnique Federale de Lausanne, Switzerland.
%
% PRACTICALSKETCHING implements a class definition for the single-view low-
% rank matrix approximation algorithms considered in [TYUC2017].
%
% PRACTICALSKETCHING is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, version 3 of the License.
%
% PRACTICALSKETCHING is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

% Find orthogonal basis "P" of the column space of "(I - UU')A"
U1tA = U1'*A;
M = A - U1*(U1tA);
P = orth(M);
% RA = P'(I - UU')A
RA = P'*M;

% Find orthogonal basis "Q" of the column space of "(I - VV')B"
V1tB = V1'*B;
N = B - V1*(V1tB);
Q = orth(N);
% RB = Q'(I - VV')B
RB = Q'*N;

% Generate matrix K in eq (4) of main reference
M = [U1tA; RA];
N = [V1tB; RB];
K = M*N';
tt = size(S1,1);
K(1:tt,1:tt) = S1 + K(1:tt,1:tt);

% Find the rotation matrices Up, Vp, and the new spectrum S
[Up, S, Vp] = svd(K,'econ');

% Generate the orthobasis U and V 
U = [U1, P]*Up;
V = [V1, Q]*Vp;

end