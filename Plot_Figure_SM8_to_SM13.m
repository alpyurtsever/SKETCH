%   This script plots the following figures in [TYUC2019]:
%   Fig.SM8, Fig.SM9, Fig.SM10, Fig.SM11, Fig.SM12, Fig.SM13
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

% Colors and more
c1 = [0,0,0];
c2 = [1,0,0];
c3 = [0,0,1];

Rsweep = [5,10,20];
for R = Rsweep

experiment = {};
experimentName = {};
isylabel = {};
islegend = {};
plotFlat = {};
xTick    = {};

% LowRank
experiment{end+1}       = ['Test_parfor_LowRankPlusNoise/HiNoise_n1000_r10_R',num2str(R),'_complex'];
experimentName{end+1}   = ['Exp2_LowRankHiNoise_R',num2str(R)];
isylabel{end+1}         = 1;
islegend{end+1}         = 0;
plotFlat{end+1}         = 1;
xTick{end+1}            = [12,24,48,96,192];

experiment{end+1}       = ['Test_parfor_LowRankPlusNoise/MedNoise_n1000_r10_R',num2str(R),'_complex'];
experimentName{end+1}   = ['Exp2_LowRankMedNoise_R',num2str(R)];
isylabel{end+1}         = 0;
islegend{end+1}         = 0;
plotFlat{end+1}         = 1;
xTick{end+1}            = xTick{end};

experiment{end+1}       = ['Test_parfor_LowRankPlusNoise/LowNoise_n1000_r10_R',num2str(R),'_complex'];
experimentName{end+1}   = ['Exp2_LowRankLowNoise_R',num2str(R)];
isylabel{end+1}         = 0;
islegend{end+1}         = 0;
plotFlat{end+1}         = 1;
xTick{end+1}            = xTick{end};

% PolyDecay
experiment{end+1}       = ['Test_parfor_PolyDecay/n1000_r10_R',num2str(R),'_p0.5_complex'];
experimentName{end+1}   = ['Exp2_PolyDecaySlow_R',num2str(R)];
isylabel{end+1}         = 1;
islegend{end+1}         = 0;
plotFlat{end+1}         = 1;
xTick{end+1}            = xTick{end};

experiment{end+1}       = ['Test_parfor_PolyDecay/n1000_r10_R',num2str(R),'_p1_complex'];
experimentName{end+1}   = ['Exp2_PolyDecayMed_R',num2str(R)];
isylabel{end+1}         = 0;
islegend{end+1}         = 0;
plotFlat{end+1}         = 1;
xTick{end+1}            = xTick{end};

experiment{end+1}       = ['Test_parfor_PolyDecay/n1000_r10_R',num2str(R),'_p2_complex'];
experimentName{end+1}   = ['Exp2_PolyDecayFast_R',num2str(R)];
isylabel{end+1}         = 0;
islegend{end+1}         = 0;
plotFlat{end+1}         = 1;
xTick{end+1}            = xTick{end};

% ExpDecay
experiment{end+1}       = ['Test_parfor_ExpDecay/n1000_r10_R',num2str(R),'_q0.01_complex'];
experimentName{end+1}   = ['Exp2_ExpDecaySlow_R',num2str(R)];
isylabel{end+1}         = 1;
islegend{end+1}         = 1;
plotFlat{end+1}         = 1;
xTick{end+1}            = xTick{end};

experiment{end+1}       = ['Test_parfor_ExpDecay/n1000_r10_R',num2str(R),'_q0.1_complex'];
experimentName{end+1}   = ['Exp2_ExpDecayMed_R',num2str(R)];
isylabel{end+1}         = 0;
islegend{end+1}         = 0;
plotFlat{end+1}         = 1;
xTick{end+1}            = xTick{end};

experiment{end+1}       = ['Test_parfor_ExpDecay/n1000_r10_R',num2str(R),'_q0.5_complex'];
experimentName{end+1}   = ['Exp2_ExpDecayFast_R',num2str(R)];
isylabel{end+1}         = 0;
islegend{end+1}         = 0;
plotFlat{end+1}         = 1;
xTick{end+1}            = xTick{end};

%% Plot S2 error
syms kk ss
assume(kk,'real');
assume(kk >= 0);
for rr = 1:length(experiment)
    clearvars Err*
    dname_gauss = dir(['results/',experiment{rr},'_Gauss/*.mat']);
    gauss = load(['results/',experiment{rr},'_Gauss/',dname_gauss.name]);
    
    numMonteCarlo = length(gauss.MCresults);
    numDataPoint  = length(gauss.MCresults{1}.info.T);
    
    r = gauss.MCresults{1}.info.r;
    RR = r;
    
    for tt = 1:numDataPoint
        T = gauss.MCresults{1}.info.T;
        
        if strcmp(gauss.MCresults{1}.info.field,'real')
            alpha = 1;
        elseif strcmp(gauss.MCresults{1}.info.field,'complex')
            alpha = 0;
        end
        
        n = gauss.MCresults{1}.info.n;
        m = gauss.MCresults{1}.info.m;
        
        k_Nat = max( r + alpha + 1, ...
            floor( (1/8) * (sqrt((m+n+4*alpha)^2 + 16*(T(tt) - alpha^2)) - (m+n+4*alpha)) ));
        
        eqns = [(T(tt) + ss^2)*(kk^2-RR^2) == 2*RR*ss^2*(ss - kk), kk*(m+n) + ss^2 == T(tt)];
        S = vpasolve(eqns, [ss,kk]);
        k_Flat = max( RR + alpha + 1, ...
            floor(S.kk)); %floor(3*R*T(tt) / (2*T(tt) + R*(m+n)));
        
        for ll = 1:numMonteCarlo
            gauss.ErrThree_S2Min(ll,tt) = min(gauss.MCresults{ll}.threesketch.ErrThree_S2{tt})     /   gauss.MCresults{ll}.ErrBest_S2     -   1;
            
            ind_Nat = find(gauss.MCresults{ll}.threesketch.k{tt} == k_Nat); % Index of theoretical (k,l) pair
            ind_Flat = find(gauss.MCresults{ll}.threesketch.k{tt} == k_Flat); % Index of theoretical (k,l) pair
            % if isempty(ind_Nat), ind_Nat = 1; end
            % if isempty(ind_Flat), ind_Flat = 1; end
            
            gauss.ErrThree_S2Nat(ll,tt)  = abs(gauss.MCresults{ll}.threesketch.ErrThree_S2{tt}(ind_Nat)    /   gauss.MCresults{ll}.ErrBest_S2     -   1); % minimum over all (k,l) pairs
            gauss.ErrThree_S2Flat(ll,tt)  = abs(gauss.MCresults{ll}.threesketch.ErrThree_S2{tt}(ind_Flat)    /   gauss.MCresults{ll}.ErrBest_S2     -   1); % minimum over all (k,l) pairs
        end
    end
    
    % Compute average relative error over all Monte Carlo iterations
    gauss.ErrThree_S2Min     = mean(gauss.ErrThree_S2Min,1);
    gauss.ErrThree_S2Nat  = mean(gauss.ErrThree_S2Nat,1);
    gauss.ErrThree_S2Flat  = mean(gauss.ErrThree_S2Flat,1);
    
    % Draw the figures now
    hfig{rr} = figure('Position',[100,100,350,320]);
    set(hfig{rr},'name',[experimentName{rr},'_S2'],'numbertitle','off')
    h = loglog(gauss.MCresults{1}.info.Tdivdim, gauss.ErrThree_S2Min,'-', 'Color',c1);
    hold on
    h = loglog(gauss.MCresults{1}.info.Tdivdim, gauss.ErrThree_S2Nat,'--', 'Color',c2);
    if plotFlat{rr}
        h = loglog(gauss.MCresults{1}.info.Tdivdim, gauss.ErrThree_S2Flat,':', 'Color',c3);
    end
    
    if islegend{rr}
        hleg = legend(...
            'Oracle', ...
            'Natural (5.6)', ...
            'Flat (5.7)' ...
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
    
    drawnow
    
end

%% Plot Sinf error
syms kk ss
assume(kk,'real');
assume(kk >= 0);
for rr = 1:length(experiment)
    clearvars Err*
    dname_gauss = dir(['results/',experiment{rr},'_Gauss/*.mat']);
    gauss = load(['results/',experiment{rr},'_Gauss/',dname_gauss.name]);
    
    numMonteCarlo = length(gauss.MCresults);
    numDataPoint  = length(gauss.MCresults{1}.info.T);
    
    r = gauss.MCresults{1}.info.r;
    RR = r;
    
    for tt = 1:numDataPoint
        T = gauss.MCresults{1}.info.T;
        
        if strcmp(gauss.MCresults{1}.info.field,'real')
            alpha = 1;
        elseif strcmp(gauss.MCresults{1}.info.field,'complex')
            alpha = 0;
        end
        
        n = gauss.MCresults{1}.info.n;
        m = gauss.MCresults{1}.info.m;
        
        k_Nat = max(r+alpha+1,...
            floor( (1/8) * (sqrt((m+n+4*alpha)^2 + 16*(T(tt) - alpha^2)) - (m+n+4*alpha)) ));
        
        eqns = [(T(tt) + ss^2)*(kk^2-RR^2) == 2*RR*ss^2*(ss - kk), kk*(m+n) + ss^2 == T(tt)];
        S = vpasolve(eqns, [ss,kk]);
        k_Flat = max(r+alpha+1,...
            floor(S.kk)); %floor(3*R*T(tt) / (2*T(tt) + R*(m+n)));
        
        for ll = 1:numMonteCarlo
            gauss.ErrThree_SinfMin(ll,tt)     = min(gauss.MCresults{ll}.threesketch.ErrThree_Sinf{tt})     /   gauss.MCresults{ll}.ErrBest_Sinf     -   1;
            
            ind_Nat = find(gauss.MCresults{ll}.threesketch.k{tt} == k_Nat); % Index of theoretical (k,l) pair
            ind_Flat = find(gauss.MCresults{ll}.threesketch.k{tt} == k_Flat); % Index of theoretical (k,l) pair
            % if isempty(ind_Nat), ind_Nat = 1; end
            % if isempty(ind_Flat), ind_Flat = 1; end
            
            gauss.ErrThree_SinfNat(ll,tt)  = abs(gauss.MCresults{ll}.threesketch.ErrThree_Sinf{tt}(ind_Nat)    /   gauss.MCresults{ll}.ErrBest_Sinf     -   1); % minimum over all (k,l) pairs
            gauss.ErrThree_SinfFlat(ll,tt)  = abs(gauss.MCresults{ll}.threesketch.ErrThree_Sinf{tt}(ind_Flat)    /   gauss.MCresults{ll}.ErrBest_Sinf     -   1); % minimum over all (k,l) pairs
        end
    end
    
    % Compute average relative error over all Monte Carlo iterations
    gauss.ErrThree_SinfMin  = mean(gauss.ErrThree_SinfMin,1);
    gauss.ErrThree_SinfNat  = mean(gauss.ErrThree_SinfNat,1);
    gauss.ErrThree_SinfFlat  = mean(gauss.ErrThree_SinfFlat,1);
    
    % Draw the figures now
    hfig{rr} = figure('Position',[100,100,350,320]);
    set(hfig{rr},'name',[experimentName{rr},'_Sinf'],'numbertitle','off')
    h = loglog(gauss.MCresults{1}.info.Tdivdim, gauss.ErrThree_SinfMin,'-', 'Color',c1);
    hold on
    h = loglog(gauss.MCresults{1}.info.Tdivdim, gauss.ErrThree_SinfNat,'--', 'Color',c2);
    if plotFlat{rr}
        h = loglog(gauss.MCresults{1}.info.Tdivdim, gauss.ErrThree_SinfFlat,':', 'Color',c3);
    end
    
    if islegend{rr}
        hleg = legend(...
            'Oracle', ...
            'Natural (5.6)', ...
            'Flat (5.7)' ...
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
    
    drawnow
    
end



end
