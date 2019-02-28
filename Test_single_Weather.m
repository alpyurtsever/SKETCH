function [sd] = Test_single_Weather(r, field, model)
%This script runs an experimental instance using the real valued Weather 
%data. See our reference paper for more details.
%	Inputs: r		-> target rank
% 			field 	-> field for the test matrices: 'complex' or 'real'.
%           model   -> model for the test matrices: 'Gauss','Sparse','SSRFT'
%	Output: sd		-> structure with following fields
%						* ErrBest_S2 ,ErrBest_S1, ErrBest_Sinf
%						* ID, twosketch, threesketch
%						* info (struct with more info about size and storage)
%						* twosketch (struct with info of rel. error by 2sketch)
%						* threesketch (struct with info of rel. error by 3sketch)
%
% 	[TYUC2019] J.A. Tropp, A. Yurtsever, M. Udell and V. Cevher.  
% 	Streaming Low-Rank Matrix Approximation with an Application to
% 	Scientific Simulation. 
%
%   Coded by: Alp Yurtsever
%   Ecole Polytechnique Federale de Lausanne, Switzerland.
%   Laboratory for Information and Inference Systems, LIONS.
%   contact: alp.yurtsever@epfl.ch
%   Last modified: June 22, 2018
%
% SKETCHv1.1
% Copyright (C) 2018 Laboratory for Information and Inference Systems
% (LIONS), Ecole Polytechnique Federale de Lausanne, Switzerland.
% This code is a part of SKETCH toolbox. 
% Please read COPYRIGHT before using this file.

%% Generate test setup and/or load data

load('data/Weather.mat');

% Find best rank-r errors
ErrBest_S2 = norm(singVals((r+1):end));
ErrBest_S1 = norm(singVals((r+1):end),1);
ErrBest_Sinf = norm(singVals((r+1):end),inf);

%% Run the generic test file
Test_generic_Type1;

end
