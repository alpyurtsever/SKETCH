%   This script plots the following figures in [TYUC2019]:
%   Fig.2, Fig.SM19
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

experiment = {};
experimentName = {};
isylabel = {};
islegend = {};
xTick    = {};

% Data
experiment{end+1}       = 'Test_parfor_Weather/r10_real_Sparse';
experimentName{end+1}   = 'Data_Weather_r10_Sparse';
isylabel{end+1}         = 1;
islegend{end+1}         = 1;
xTick{end+1}            = 12.*(2.^(0:10)');
 
experiment{end+1}       = 'Test_parfor_NavierStokes/r10_StreamVelocity_real_Sparse';
experimentName{end+1}   = 'Data_DNS_StreamVelocity_r10_Sparse';
isylabel{end+1}         = 0;
islegend{end+1}         = 0;
xTick{end+1}            = xTick{end};

experiment{end+1}       = 'Test_parfor_MaxCut/r1_real_Sparse';
experimentName{end+1}   = 'Data_MaxCut_r1_Sparse';
isylabel{end+1}         = 1;
islegend{end+1}         = 0;
xTick{end+1}            = xTick{end};

experiment{end+1}       = 'Test_parfor_MaxCut/r14_real_Sparse';
experimentName{end+1}   = 'Data_MaxCut_r14_Sparse';
isylabel{end+1}         = 0;
islegend{end+1}         = 0;
xTick{end+1}            = xTick{end};

experiment{end+1}       = 'Test_parfor_PhaseRetrieval/r1_complex_Sparse';
experimentName{end+1}   = 'Data_PhaseRetrieval_r1_Sparse';
isylabel{end+1}         = 1;
islegend{end+1}         = 0;
xTick{end+1}            = 12.*(2.^(0:10)');

experiment{end+1}       = 'Test_parfor_PhaseRetrieval/r5_complex_Sparse';
experimentName{end+1}   = 'Data_PhaseRetrieval_r5_Sparse';
isylabel{end+1}         = 0;
islegend{end+1}         = 0;
xTick{end+1}            = xTick{end};

% Choose colors
cAlgBetter  = [0,0,0];
cAlgBetter2 = [0.6,0.6,0.6];
cAlg7    = [.6,.1,0];
cAlg72   = [.75,.4,.35];
cAlgUpa  = [0,0,.9];
cAlgUpa2 = [0,.7,.9];
cAlgHmt  = [ 0.900 0.760 0.0700];
cAlgHmt2  = [1,0.9,0.25];

%% Plot S2 error
for rr = 1:length(experiment)
    clearvars Err*
    dname_gauss = dir(['results/',experiment{rr},'/*.mat']);
    gauss = load(['results/',experiment{rr},'/',dname_gauss.name]);

    numMonteCarlo = length(gauss.MCresults);
    numDataPoint  = length(gauss.MCresults{1}.info.T);
    
    for ll = 1:numMonteCarlo
        for tt = 1:numDataPoint
            gauss.ErrThree_S2Min(ll,tt)     = min(gauss.MCresults{ll}.threesketch.ErrThree_S2{tt})     /   gauss.MCresults{ll}.ErrBest_S2     -   1; 
            if isempty(gauss.MCresults{ll}.twosketch.ErrA7_S2{tt}), gauss.MCresults{ll}.twosketch.ErrA7_S2{tt} = nan; end
            gauss.ErrA7_S2Min(ll,tt)     = min(gauss.MCresults{ll}.twosketch.ErrA7_S2{tt})     /   gauss.MCresults{ll}.ErrBest_S2     -   1; % minimum over all (k,l) pairs
            if isempty(gauss.MCresults{ll}.threesketch.ErrUpa_S2{tt}), gauss.MCresults{ll}.threesketch.ErrUpa_S2{tt} = nan; end
            gauss.ErrUpa_S2Min(ll,tt)     = min(gauss.MCresults{ll}.threesketch.ErrUpa_S2{tt})     /   gauss.MCresults{ll}.ErrBest_S2     -   1; 
            gauss.ErrHmt_S2Min(ll,tt) = gauss.MCresults{ll}.hmtsketch.ErrHMT_S2(tt) / gauss.MCresults{ll}.ErrBest_S2  -   1; 
            
            T = gauss.MCresults{ll}.info.T;
            if strcmp(gauss.MCresults{ll}.info.field,'real')
                alpha = 1;
            elseif strcmp(gauss.MCresults{ll}.info.field,'complex')
                alpha = 0;
            else
                error('field should be ''real'' or ''complex''');
            end
            n = gauss.MCresults{ll}.info.n;
            m = gauss.MCresults{ll}.info.m;

            k_Nat = floor( (1/8) * (sqrt((m+n+4*alpha)^2 + 16*(T(tt) - alpha^2)) - (m+n+4*alpha)) );
            ind_Nat = find(gauss.MCresults{ll}.threesketch.k{tt} == k_Nat); % Index of theoretical (k,l) pair
            if isempty(ind_Nat), ind_Nat = 1; end
            gauss.ErrThree_S2Nat(ll,tt)  = abs(gauss.MCresults{ll}.threesketch.ErrThree_S2{tt}(ind_Nat) / gauss.MCresults{ll}.ErrBest_S2     -   1); % minimum over all (k,l) pairs
            gauss.ErrUpa_S2Nat(ll,tt)  = abs(gauss.MCresults{ll}.threesketch.ErrUpa_S2{tt}(ind_Nat) / gauss.MCresults{ll}.ErrBest_S2     -   1); % minimum over all (k,l) pairs
        
            Tdivdim = gauss.MCresults{ll}.info.Tdivdim(tt);
            r = gauss.MCresults{ll}.info.r;
            
            k_Decay_A7 = max(r+alpha+1,floor((T(tt) - n*alpha)/(m+2*n)));

            ind_Decay_A7 = find(gauss.MCresults{ll}.twosketch.k{tt} == k_Decay_A7); % Index of theoretical (k,l) pair
            if isempty(ind_Decay_A7), ind_Decay_A7 = 1; end
            gauss.ErrA7_S2Decay(ll,tt)  = abs(gauss.MCresults{ll}.twosketch.ErrA7_S2{tt}(ind_Decay_A7) / gauss.MCresults{ll}.ErrBest_S2     -   1); % minimum over all (k,l) pairs
        end
    end
    
    % Compute average relative error over all Monte Carlo iterations
    gauss.ErrThree_S2Min     = mean(gauss.ErrThree_S2Min,1);
    gauss.ErrA7_S2Min     = mean(gauss.ErrA7_S2Min,1);
    gauss.ErrUpa_S2Min     = mean(gauss.ErrUpa_S2Min,1);
    gauss.ErrThree_S2Nat     = mean(gauss.ErrThree_S2Nat,1);
    gauss.ErrA7_S2Decay     = mean(gauss.ErrA7_S2Decay,1);
    gauss.ErrUpa_S2Nat     = mean(gauss.ErrUpa_S2Nat,1);
    gauss.ErrHmt_S2Min = mean(gauss.ErrHmt_S2Min,1);
    
    % Draw the figures now
    hfig{rr} = figure('Position',[100,100,350,320]);
    set(hfig{rr},'name',[experimentName{rr},'_S2'],'numbertitle','off')
    h1 = loglog(gauss.MCresults{1}.info.Tdivdim, gauss.ErrA7_S2Min,'-', 'Color',cAlg7);
    hold on
    h2 = loglog(gauss.MCresults{1}.info.Tdivdim, gauss.ErrUpa_S2Min,'-', 'Color',cAlgUpa);
    h3 = loglog(gauss.MCresults{1}.hmtsketch.l, gauss.ErrHmt_S2Min,'-', 'Color',cAlgHmt);
    h4 = loglog(gauss.MCresults{1}.info.Tdivdim, gauss.ErrThree_S2Min,'-', 'Color',cAlgBetter);
    h = loglog(gauss.MCresults{1}.info.Tdivdim, gauss.ErrA7_S2Decay,'--', 'Color',cAlg72);
    h = loglog(gauss.MCresults{1}.info.Tdivdim, gauss.ErrUpa_S2Nat,'--', 'Color',cAlgUpa2);
    h = loglog(gauss.MCresults{1}.info.Tdivdim, gauss.ErrThree_S2Nat,'--', 'Color',cAlgBetter2);
    
    % Add the legend
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
    
%% Plot Sinf error

for rr = 1:length(experiment)
    clearvars Err*
    dname_gauss = dir(['results/',experiment{rr},'/*.mat']);
    gauss = load(['results/',experiment{rr},'/',dname_gauss.name]);

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
            
            T = gauss.MCresults{ll}.info.T;
            if strcmp(gauss.MCresults{ll}.info.field,'real')
                alpha = 1;
            elseif strcmp(gauss.MCresults{ll}.info.field,'complex')
                alpha = 0;
            else
                error('field should be ''real'' or ''complex''');
            end
            n = gauss.MCresults{ll}.info.n;
            m = gauss.MCresults{ll}.info.m;

            k_Nat = floor( (1/8) * (sqrt((m+n+4*alpha)^2 + 16*(T(tt) - alpha^2)) - (m+n+4*alpha)) );
            ind_Nat = find(gauss.MCresults{ll}.threesketch.k{tt} == k_Nat); % Index of theoretical (k,l) pair
            if isempty(ind_Nat), ind_Nat = 1; end
            gauss.ErrThree_SinfNat(ll,tt)  = abs(gauss.MCresults{ll}.threesketch.ErrThree_Sinf{tt}(ind_Nat) / gauss.MCresults{ll}.ErrBest_Sinf     -   1); % minimum over all (k,l) pairs
            gauss.ErrUpa_SinfNat(ll,tt)  = abs(gauss.MCresults{ll}.threesketch.ErrUpa_Sinf{tt}(ind_Nat) / gauss.MCresults{ll}.ErrBest_Sinf     -   1); % minimum over all (k,l) pairs
        
            Tdivdim = gauss.MCresults{ll}.info.Tdivdim(tt);
            r = gauss.MCresults{ll}.info.r;
            
            k_Decay_A7 = max(r+alpha+1,floor((T(tt) - n*alpha)/(m+2*n)));

            ind_Decay_A7 = find(gauss.MCresults{ll}.twosketch.k{tt} == k_Decay_A7); % Index of theoretical (k,l) pair
            if isempty(ind_Decay_A7), ind_Decay_A7 = 1; end
            gauss.ErrA7_SinfDecay(ll,tt)  = abs(gauss.MCresults{ll}.twosketch.ErrA7_Sinf{tt}(ind_Decay_A7) / gauss.MCresults{ll}.ErrBest_Sinf     -   1); % minimum over all (k,l) pairs
        end
    end
    
    % Compute average relative error over all Monte Carlo iterations
    gauss.ErrThree_SinfMin     = mean(gauss.ErrThree_SinfMin,1);
    gauss.ErrA7_SinfMin     = mean(gauss.ErrA7_SinfMin,1);
    gauss.ErrUpa_SinfMin     = mean(gauss.ErrUpa_SinfMin,1);
    gauss.ErrThree_SinfNat     = mean(gauss.ErrThree_SinfNat,1);
    gauss.ErrA7_SinfDecay     = mean(gauss.ErrA7_SinfDecay,1);
    gauss.ErrUpa_SinfNat     = mean(gauss.ErrUpa_SinfNat,1);
    gauss.ErrHmt_SinfMin = mean(gauss.ErrHmt_SinfMin,1);

    % Draw the figures now
    hfig{rr} = figure('Position',[100,100,350,320]);
    set(hfig{rr},'name',[experimentName{rr},'_Sinf'],'numbertitle','off')
    h1 = loglog(gauss.MCresults{1}.info.Tdivdim, gauss.ErrA7_SinfMin,'-', 'Color',cAlg7);
    hold on
    h2 = loglog(gauss.MCresults{1}.info.Tdivdim, gauss.ErrUpa_SinfMin,'-', 'Color',cAlgUpa);
    h3 = loglog(gauss.MCresults{1}.hmtsketch.l, gauss.ErrHmt_SinfMin,'-', 'Color',cAlgHmt);
    h4 = loglog(gauss.MCresults{1}.info.Tdivdim, gauss.ErrThree_SinfMin,'-', 'Color',cAlgBetter);
    h = loglog(gauss.MCresults{1}.info.Tdivdim, gauss.ErrA7_SinfDecay,'--', 'Color',cAlg72);
    h = loglog(gauss.MCresults{1}.info.Tdivdim, gauss.ErrUpa_SinfNat,'--', 'Color',cAlgUpa2);
    h = loglog(gauss.MCresults{1}.info.Tdivdim, gauss.ErrThree_SinfNat,'--', 'Color',cAlgBetter2);
    
    % Add the legend
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
