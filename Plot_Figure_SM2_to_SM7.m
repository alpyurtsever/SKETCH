%   This script plots the following figures in [TYUC2019]:
%   Fig.SM2, Fig.SM3, Fig.SM4, Fig.SM5, Fig.SM6, Fig.SM7
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

% Colors 
c1 = [0,0,0];
c2 = [0,0,1];
c3 = [1,0,0];

experiment = {};
experimentName = {};
isylabel = {};
islegend = {};
plotFlat = {};
xTick    = {};

Rsweep = [5,10,20];
for R = Rsweep
% Data
experiment{end+1}       = ['Test_parfor_LowRankPlusNoise/HiNoise_n1000_r10_R',num2str(R),'_complex'];
experimentName{end+1}   = ['Exp3_LowRankHiNoise_R',num2str(R)];
isylabel{end+1}         = 1;
islegend{end+1}         = 0;
plotFlat{end+1}         = 1;
xTick{end+1}            = [12,24,48,96,192];

experiment{end+1}       = ['Test_parfor_LowRankPlusNoise/MedNoise_n1000_r10_R',num2str(R),'_complex'];
experimentName{end+1}   = ['Exp3_LowRankMedNoise_R',num2str(R)];
isylabel{end+1}         = 0;
islegend{end+1}         = 0;
plotFlat{end+1}         = 1;
xTick{end+1}            = xTick{end};

experiment{end+1}       = ['Test_parfor_LowRankPlusNoise/LowNoise_n1000_r10_R',num2str(R),'_complex'];
experimentName{end+1}   = ['Exp3_LowRankLowNoise_R',num2str(R)];
isylabel{end+1}         = 0;
islegend{end+1}         = 0;
plotFlat{end+1}         = 1;
xTick{end+1}            = xTick{end};

% Low-Rank + poly decay

experiment{end+1}       = ['Test_parfor_PolyDecay/n1000_r10_R',num2str(R),'_p0.5_complex'];
experimentName{end+1}   = ['Exp3_PolyDecaySlow_R',num2str(R)];
isylabel{end+1}         = 1;
islegend{end+1}         = 0;
plotFlat{end+1}         = 0;
xTick{end+1}            = xTick{end};


experiment{end+1}       = ['Test_parfor_PolyDecay/n1000_r10_R',num2str(R),'_p1_complex'];
experimentName{end+1}   = ['Exp3_PolyDecayMed_R',num2str(R)];
isylabel{end+1}         = 0;
islegend{end+1}         = 0;
plotFlat{end+1}         = 0;
xTick{end+1}            = xTick{end};

experiment{end+1}       = ['Test_parfor_PolyDecay/n1000_r10_R',num2str(R),'_p2_complex'];
experimentName{end+1}   = ['Exp3_PolyDecayFast_R',num2str(R)];
isylabel{end+1}         = 0;
islegend{end+1}         = 0;
plotFlat{end+1}         = 0;
xTick{end+1}            = xTick{end};

% Low-Rank + exp decay

experiment{end+1}       = ['Test_parfor_ExpDecay/n1000_r10_R',num2str(R),'_q0.01_complex'];
experimentName{end+1}   = ['Exp3_ExpDecaySlow_R',num2str(R)];
isylabel{end+1}         = 1;
islegend{end+1}         = 1;
plotFlat{end+1}         = 0;
xTick{end+1}            = xTick{end};

experiment{end+1}       = ['Test_parfor_ExpDecay/n1000_r10_R',num2str(R),'_q0.1_complex'];
experimentName{end+1}   = ['Exp3_ExpDecayMed_R',num2str(R)];
isylabel{end+1}         = 0;
islegend{end+1}         = 0;
plotFlat{end+1}         = 0;
xTick{end+1}            = xTick{end};

experiment{end+1}       =['Test_parfor_ExpDecay/n1000_r10_R',num2str(R),'_q0.5_complex'];
experimentName{end+1}   = ['Exp3_ExpDecayFast_R',num2str(R)];
isylabel{end+1}         = 0;
islegend{end+1}         = 0;
plotFlat{end+1}         = 0;
xTick{end+1}            = xTick{end};

end

%% Plot S2 error
for rr = 1:length(experiment)
    clearvars Err*
    dname_gauss = dir(['results/',experiment{rr},'_Gauss/*.mat']);
    gauss = load(['results/',experiment{rr},'_Gauss/',dname_gauss.name]);

    numMonteCarlo = length(gauss.MCresults);
    numDataPoint  = length(gauss.MCresults{1}.info.T);
    
    for ll = 1:numMonteCarlo
        for tt = 1:numDataPoint
            gauss.ErrThree_S2Min(ll,tt)     = min(gauss.MCresults{ll}.threesketch.ErrThree_S2{tt})     /   gauss.MCresults{ll}.ErrBest_S2     -   1; 
        end
    end
    
    % Compute average relative error over all Monte Carlo iterations
    gauss.ErrThree_S2Min     = mean(gauss.ErrThree_S2Min,1);
    
    dname_ssrft = dir(['results/',experiment{rr},'_SSRFT/*.mat']);
    ssrft = load(['results/',experiment{rr},'_SSRFT/',dname_ssrft.name]);

    numMonteCarlo = length(ssrft.MCresults);
    numDataPoint  = length(ssrft.MCresults{1}.info.T);
    
    for ll = 1:numMonteCarlo
        for tt = 1:numDataPoint
            ssrft.ErrThree_S2Min(ll,tt)     = min(ssrft.MCresults{ll}.threesketch.ErrThree_S2{tt})     /   gauss.MCresults{ll}.ErrBest_S2     -   1; 
        end
    end
    
    % Compute average relative error over all Monte Carlo iterations
    ssrft.ErrThree_S2Min     = mean(ssrft.ErrThree_S2Min,1);
    
    dname_sparse = dir(['results/',experiment{rr},'_Sparse/*.mat']);
    spars = load(['results/',experiment{rr},'_Sparse/',dname_sparse.name]);

    numMonteCarlo = length(spars.MCresults);
    numDataPoint  = length(spars.MCresults{1}.info.T);
    
    for ll = 1:numMonteCarlo
        for tt = 1:numDataPoint
            spars.ErrThree_S2Min(ll,tt)     = min(spars.MCresults{ll}.threesketch.ErrThree_S2{tt})     /   gauss.MCresults{ll}.ErrBest_S2     -   1; 
        end
    end
    
    % Compute average relative error over all Monte Carlo iterations
    spars.ErrThree_S2Min     = mean(spars.ErrThree_S2Min,1);
    
    
    % Draw the figures now
    hfig{rr} = figure('Position',[100,100,350,320]);
    set(hfig{rr},'name',[experimentName{rr},'_S2'],'numbertitle','off')
    h = loglog(gauss.MCresults{1}.info.Tdivdim, gauss.ErrThree_S2Min,'-', 'Color',c1);
    hold on
    h = loglog(ssrft.MCresults{1}.info.Tdivdim, ssrft.ErrThree_S2Min,'--', 'Color',c2);
    h = loglog(spars.MCresults{1}.info.Tdivdim, spars.ErrThree_S2Min,':', 'Color',c3);

    if islegend{rr}
    hleg = legend(...
        'Gauss', ...
        'SSRFT', ...
        'Sparse' ...
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
    xlim([gauss.MCresults{1}.info.Tdivdim(1),gauss.MCresults{1}.info.Tdivdim(end)])
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
        end
    end
    
    % Compute average relative error over all Monte Carlo iterations
    gauss.ErrThree_SinfMin     = mean(gauss.ErrThree_SinfMin,1);
    
    dname_ssrft = dir(['results/',experiment{rr},'_SSRFT/*.mat']);
    ssrft = load(['results/',experiment{rr},'_SSRFT/',dname_ssrft.name]);

    numMonteCarlo = length(ssrft.MCresults);
    numDataPoint  = length(ssrft.MCresults{1}.info.T);
    
    for ll = 1:numMonteCarlo
        for tt = 1:numDataPoint
            ssrft.ErrThree_SinfMin(ll,tt)     = min(ssrft.MCresults{ll}.threesketch.ErrThree_Sinf{tt})     /   gauss.MCresults{ll}.ErrBest_Sinf     -   1; 
        end
    end
    
    % Compute average relative error over all Monte Carlo iterations
    ssrft.ErrThree_SinfMin     = mean(ssrft.ErrThree_SinfMin,1);
    
    dname_sparse = dir(['results/',experiment{rr},'_Sparse/*.mat']);
    spars = load(['results/',experiment{rr},'_Sparse/',dname_sparse.name]);

    numMonteCarlo = length(spars.MCresults);
    numDataPoint  = length(spars.MCresults{1}.info.T);
    
    for ll = 1:numMonteCarlo
        for tt = 1:numDataPoint
            spars.ErrThree_SinfMin(ll,tt)     = min(spars.MCresults{ll}.threesketch.ErrThree_Sinf{tt})     /   gauss.MCresults{ll}.ErrBest_Sinf     -   1; 
        end
    end
    
    % Compute average relative error over all Monte Carlo iterations
    spars.ErrThree_SinfMin     = mean(spars.ErrThree_SinfMin,1);
    
    
    % Draw the figures now
    hfig{rr} = figure('Position',[100,100,350,320]);
    set(hfig{rr},'name',[experimentName{rr},'_Sinf'],'numbertitle','off')
    h = loglog(gauss.MCresults{1}.info.Tdivdim, gauss.ErrThree_SinfMin,'-', 'Color',c1);
    hold on
    h = loglog(ssrft.MCresults{1}.info.Tdivdim, ssrft.ErrThree_SinfMin,'--', 'Color',c2);
    h = loglog(spars.MCresults{1}.info.Tdivdim, spars.ErrThree_SinfMin,':', 'Color',c3);

    if islegend{rr}
    hleg = legend(...
        'Gauss', ...
        'SSRFT', ...
        'Sparse' ...
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
    xlim([gauss.MCresults{1}.info.Tdivdim(1),gauss.MCresults{1}.info.Tdivdim(end)])
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


