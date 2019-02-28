%   This script plots the following figures in [TYUC2019]:
%   Fig.1, Fig.SM14, Fig.SM15, Fig.SM16, Fig.SM17, Fig.SM18
%
%	[TYUC2019] J.A. Tropp, A. Yurtsever, M. Udell and V. Cevher.  
%	Streaming Low-Rank Matrix Approximation with an Application to
%	Scientific Simulation. 
%
%   Coded by: Alp Yurtsever
%   Ecole Polytechnique Federale de Lausanne, Switzerland.
%   Laboratory for Information and Inference Systems, LIONS.
%   contact: alp.yurtsever@epfl.ch
%   Last modified: February 22, 2019
%
% SKETCHv1.1
% Copyright (C) 2018 Laboratory for Information and Inference Systems
% (LIONS), Ecole Polytechnique Federale de Lausanne, Switzerland.
% This code is a part of SKETCH toolbox. 
% Please read COPYRIGHT before using this file.

%% LOAD DATA
clearvars
close all
hfig = {};

% Choose colors
cAlgBetter  = [0,0,0];
cAlg7    = [.6,.1,0];
cAlgUpa  = [0,0,.9];
cAlgWlrt = [0.9,0.6,0.15];
cAlgHmt  = [ 0.900 0.760 0.0700];

% Values of R
Rsweep = [5,10,20];
for R = Rsweep
    
experiment = {};
experimentName = {};
isylabel = {};
islegend = {};
xTick    = {};

% LowRank
experiment{end+1}       = ['Test_parfor_LowRankPlusNoise/HiNoise_n1000_r10_R',num2str(R),'_complex'];
experimentName{end+1}   = ['Exp1_LowRankHiNoise_R',num2str(R)];
isylabel{end+1}         = 1;
islegend{end+1}         = 0;
xTick{end+1}            = [12,24,48,96,192];

experiment{end+1}       = ['Test_parfor_LowRankPlusNoise/MedNoise_n1000_r10_R',num2str(R),'_complex'];
experimentName{end+1}   = ['Exp1_LowRankMedNoise_R',num2str(R)];
isylabel{end+1}         = 0;
islegend{end+1}         = 0;
xTick{end+1}            = xTick{end};

experiment{end+1}       = ['Test_parfor_LowRankPlusNoise/LowNoise_n1000_r10_R',num2str(R),'_complex'];
experimentName{end+1}   = ['Exp1_LowRankLowNoise_R',num2str(R)];
isylabel{end+1}         = 0;
islegend{end+1}         = 0;
xTick{end+1}            = xTick{end};

% PolyDecay
experiment{end+1}       = ['Test_parfor_PolyDecay/n1000_r10_R',num2str(R),'_p0.5_complex'];
experimentName{end+1}   = ['Exp1_PolyDecaySlow_R',num2str(R)];
isylabel{end+1}         = 1;
islegend{end+1}         = 0;
xTick{end+1}            = xTick{end};

experiment{end+1}       = ['Test_parfor_PolyDecay/n1000_r10_R',num2str(R),'_p1_complex'];
experimentName{end+1}   = ['Exp1_PolyDecayMed_R',num2str(R)];
isylabel{end+1}         = 0;
islegend{end+1}         = 0;
xTick{end+1}            = xTick{end};

experiment{end+1}       = ['Test_parfor_PolyDecay/n1000_r10_R',num2str(R),'_p2_complex'];
experimentName{end+1}   = ['Exp1_PolyDecayFast_R',num2str(R)];
isylabel{end+1}         = 0;
islegend{end+1}         = 0;
xTick{end+1}            = xTick{end};

% ExpDecay
experiment{end+1}       = ['Test_parfor_ExpDecay/n1000_r10_R',num2str(R),'_q0.01_complex'];
experimentName{end+1}   = ['Exp1_ExpDecaySlow_R',num2str(R)];
isylabel{end+1}         = 1;
islegend{end+1}         = 1;
xTick{end+1}            = xTick{end};

experiment{end+1}       = ['Test_parfor_ExpDecay/n1000_r10_R',num2str(R),'_q0.1_complex'];
experimentName{end+1}   = ['Exp1_ExpDecayMed_R',num2str(R)];
isylabel{end+1}         = 0;
islegend{end+1}         = 0;
xTick{end+1}            = xTick{end};

experiment{end+1}       = ['Test_parfor_ExpDecay/n1000_r10_R',num2str(R),'_q0.5_complex'];
experimentName{end+1}   = ['Exp1_ExpDecayFast_R',num2str(R)];
isylabel{end+1}         = 0;
islegend{end+1}         = 0;
xTick{end+1}            = xTick{end};

%% Plot S2 error
for rr = 1:length(experiment)
    clearvars Err*
    dname_gauss = dir(['results/',experiment{rr},'_Gauss/*.mat']);
    gauss = load(['results/',experiment{rr},'_Gauss/',dname_gauss.name]);

    numMonteCarlo = length(gauss.MCresults);
    numDataPoint  = length(gauss.MCresults{1}.info.T);
    
    for ll = 1:numMonteCarlo
        for tt = 1:numDataPoint
            gauss.ErrThree_S2Min(ll,tt)  = min(gauss.MCresults{ll}.threesketch.ErrThree_S2{tt})     /   gauss.MCresults{ll}.ErrBest_S2     -   1; 
            if isempty(gauss.MCresults{ll}.twosketch.ErrA7_S2{tt}), gauss.MCresults{ll}.twosketch.ErrA7_S2{tt} = nan; end
            gauss.ErrA7_S2Min(ll,tt)     = min(gauss.MCresults{ll}.twosketch.ErrA7_S2{tt})     /   gauss.MCresults{ll}.ErrBest_S2     -   1; % minimum over all (k,l) pairs
            if isempty(gauss.MCresults{ll}.threesketch.ErrUpa_S2{tt}), gauss.MCresults{ll}.threesketch.ErrUpa_S2{tt} = nan; end
            gauss.ErrUpa_S2Min(ll,tt)     = min(gauss.MCresults{ll}.threesketch.ErrUpa_S2{tt})     /   gauss.MCresults{ll}.ErrBest_S2     -   1; 
            gauss.ErrHmt_S2Min(ll,tt) = gauss.MCresults{ll}.hmtsketch.ErrHMT_S2(tt) / gauss.MCresults{ll}.ErrBest_S2  -   1; 
        end
    end
    [gauss.MCresults{1}.hmtsketch.l, ic] = unique(gauss.MCresults{1}.hmtsketch.l);
    gauss.ErrHmt_S2Min = gauss.ErrHmt_S2Min(:,ic);
    
    % Compute average relative error over all Monte Carlo iterations
    gauss.ErrThree_S2Min = mean(gauss.ErrThree_S2Min,1);
    gauss.ErrA7_S2Min = mean(gauss.ErrA7_S2Min,1);
    gauss.ErrUpa_S2Min = mean(gauss.ErrUpa_S2Min,1);
    gauss.ErrHmt_S2Min = mean(gauss.ErrHmt_S2Min,1);
    
    % Draw the figures now
    hfig{rr} = figure('Position',[100,100,350,320]);
    set(hfig{rr},'name',[experimentName{rr},'_S2'],'numbertitle','off')
    h1 = loglog(gauss.MCresults{1}.info.Tdivdim, gauss.ErrA7_S2Min,'-', 'Color',cAlg7);
    hold on
    h2 = loglog(gauss.MCresults{1}.info.Tdivdim, gauss.ErrUpa_S2Min,'-', 'Color',cAlgUpa);
    h3 = loglog(gauss.MCresults{1}.hmtsketch.l, gauss.ErrHmt_S2Min,'-', 'Color',cAlgHmt);
    h4 = loglog(gauss.MCresults{1}.info.Tdivdim, gauss.ErrThree_S2Min,'-', 'Color',cAlgBetter);

    if islegend{rr}
    hleg = legend([h3,h2,h1,h4],...
        '[HMT11]', ...
        '[UPA16]', ...
        '[TYUC17]', ...
        'Eqn. $(2.10)$' ...
        );
    hleg.FontSize = 16;
    hleg.Interpreter = 'latex';
    hleg.Location = 'SouthWest';
    end
    
    ax = gca;
    ax.XTick = xTick{rr};
    set(ax, 'FontSize', 13)
    xlabel('Storage:  $T/(m+n)$','Interpreter','latex','FontSize',17);
    if isylabel{rr}
    ylabel('Relative Error ($S_2$)','Interpreter','latex','FontSize',19);
    end
    
    axis tight
    xlim([-inf, inf])
    ylimCur = ylim;
    ylimCur(1) = max(ylimCur(1),1e-9);
    if ylimCur(2)/100 < ylimCur(1)
        ylimCur(1) = 0.9*10^(floor(log10(ylimCur(1))));
    end
    ylim(ylimCur);
    set(gca,'TickDir','out')
    set(gca,'LineWidth',1,'TickLength',[0.02 0.02]);
    set(findall(gca, 'Type', 'Line'),'LineWidth',2.5); 

    if rr == 1
        setPos = ax.Position;
    else
        ax.Position = setPos;
    end
    
end

%% Plot Sinf error
for rr = 1:length(experiment)
    clearvars Err*
    dname_gauss = dir(['results/',experiment{rr},'_Gauss/*.mat']);
    gauss = load(['results/',experiment{rr},'_Gauss/',dname_gauss.name]);

    numMonteCarlo = length(gauss.MCresults);
    numDataPoint  = length(gauss.MCresults{1}.info.T);
    
    for ll = 1:numMonteCarlo
        for tt = 1:numDataPoint
            gauss.ErrThree_SinfMin(ll,tt)     = min(gauss.MCresults{ll}.threesketch.ErrThree_Sinf{tt})     /   gauss.MCresults{ll}.ErrBest_Sinf     -   1; 
            if isempty(gauss.MCresults{ll}.twosketch.ErrA7_Sinf{tt}), gauss.MCresults{ll}.twosketch.ErrA7_Sinf{tt} = nan; end
            gauss.ErrA7_SinfMin(ll,tt)     = min(gauss.MCresults{ll}.twosketch.ErrA7_Sinf{tt})     /   gauss.MCresults{ll}.ErrBest_Sinf     -   1; % minimum over all (k,l) pairs
            if isempty(gauss.MCresults{ll}.threesketch.ErrUpa_Sinf{tt}), gauss.MCresults{ll}.threesketch.ErrUpa_Sinf{tt} = nan; end
            gauss.ErrUpa_SinfMin(ll,tt)     = min(gauss.MCresults{ll}.threesketch.ErrUpa_Sinf{tt})     /   gauss.MCresults{ll}.ErrBest_Sinf     -   1; 
            gauss.ErrHmt_SinfMin(ll,tt) = gauss.MCresults{ll}.hmtsketch.ErrHMT_Sinf(tt) / gauss.MCresults{ll}.ErrBest_Sinf  -   1; 
        end
    end
    [gauss.MCresults{1}.hmtsketch.l, ic] = unique(gauss.MCresults{1}.hmtsketch.l);
    gauss.ErrHmt_SinfMin = gauss.ErrHmt_SinfMin(:,ic);

    % Compute average relative error over all Monte Carlo iterations
    gauss.ErrThree_SinfMin = mean(gauss.ErrThree_SinfMin,1);
    gauss.ErrA7_SinfMin = mean(gauss.ErrA7_SinfMin,1);
    gauss.ErrUpa_SinfMin = mean(gauss.ErrUpa_SinfMin,1);
    gauss.ErrHmt_SinfMin = mean(gauss.ErrHmt_SinfMin,1);

    % Draw the figures now
    hfig{rr} = figure('Position',[100,100,350,320]);
    set(hfig{rr},'name',[experimentName{rr},'_Sinf'],'numbertitle','off')
    h1 = loglog(gauss.MCresults{1}.info.Tdivdim, gauss.ErrA7_SinfMin,'-', 'Color',cAlg7);
    hold on
    h2 = loglog(gauss.MCresults{1}.info.Tdivdim, gauss.ErrUpa_SinfMin,'-', 'Color',cAlgUpa);
    h3 = loglog(gauss.MCresults{1}.hmtsketch.l, gauss.ErrHmt_SinfMin,'-', 'Color',cAlgHmt);
    h4 = loglog(gauss.MCresults{1}.info.Tdivdim, gauss.ErrThree_SinfMin,'-', 'Color',cAlgBetter);
    
    if islegend{rr}
    hleg = legend([h3,h2,h1,h4],...
        '[HMT11]', ...
        '[UPA16]', ...
        '[TYUC17]', ...
        'Eqn. $(2.10)$' ...
        );
    hleg.FontSize = 16;
    hleg.Interpreter = 'latex';
    hleg.Location = 'SouthWest';
    end
    
    ax = gca;
    ax.XTick = xTick{rr};
    set(ax, 'FontSize', 13)
    xlabel('Storage:  $T/(m+n)$','Interpreter','latex','FontSize',17);
    if isylabel{rr}
    ylabel('Relative Error ($S_{\infty}$)','Interpreter','latex','FontSize',19);
    end
    
    axis tight
    xlim([-inf, inf])
    ylimCur = ylim;
    ylimCur(1) = max(ylimCur(1),1e-9);
    if ylimCur(2)/100 < ylimCur(1)
        ylimCur(1) = 0.9*10^(floor(log10(ylimCur(1))));
    end
    ylim(ylimCur);
    set(gca,'TickDir','out')
    set(gca,'LineWidth',1,'TickLength',[0.02 0.02]);
    set(findall(gca, 'Type', 'Line'),'LineWidth',2.5); 

    if rr == 1
        setPos = ax.Position;
    else
        ax.Position = setPos;
    end
    
end

end
