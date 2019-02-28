function [dummy] = Test_parfor_MaxCut(r,field,model,MonteCarloNum)
%This function runs Monte Carlo simulations in a parfor loop.
%
% [TYUC2019] J.A. Tropp, A. Yurtsever, M. Udell and V. Cevher.  
% Streaming Low-Rank Matrix Approximation with an Application to
% Scientific Simulation. 
%
% Coded by: Alp Yurtsever
% Ecole Polytechnique Federale de Lausanne, Switzerland.
% Laboratory for Information and Inference Systems, LIONS.
% contact: alp.yurtsever@epfl.ch
% Last modified: June 22, 2018
%
% SKETCHv1.1
% Copyright (C) 2018 Laboratory for Information and Inference Systems
% (LIONS), Ecole Polytechnique Federale de Lausanne, Switzerland.
% This code is a part of SKETCH toolbox. 
% Please read COPYRIGHT before using this file.

MCresults = cell(MonteCarloNum,1);

parfor tt = 1:MonteCarloNum
    MCresults{tt} = Test_single_MaxCut(r,field,model);
end

name = ['r',num2str(r),'_',field,'_',model,'/'];
savepath = ['results/',mfilename,'/',name];
if ~exist(savepath,'dir'), mkdir(savepath); end
filename = datestr(now,30);
save([savepath,filename],'-v7.3','MCresults');

dummy = 1;

end