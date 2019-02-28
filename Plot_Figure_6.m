%   This script plots the following figures in [TYUC2019]:
%   Fig.6
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

%% Colors and more
close all
dname = dir('results/Test_parfor_ErrorEstimate_NavierStokes/*.mat');
load(['results/Test_parfor_ErrorEstimate_NavierStokes/',dname.name]);

%% Figure 1
%
% close all

ErrBest = sqrt(cumsum(MCresults{1}.spectrumExact(1:end).^2,'reverse'));
ErrBest = circshift(ErrBest,[-1,0]);
ErrBest(end) = 0;

%% Plot for k = 16
k=16;
hfig{1} = figure('Position',[100,100,350,320]);
set(hfig{1},'name','scatter-approximation-k=16','numbertitle','off');
xx = nan(length(MCresults),1);
yy = nan(length(MCresults),1);
for t = 1:length(MCresults)
    xx(t) = MCresults{t}.ErrLowRank(k).^2;
    yy(t) = MCresults{t}.ErrLowRankEstimate{2}(k).^2;
end
hold on
hs = scatter(xx,yy);
hs.MarkerFaceColor = [0,0,1];
hs.MarkerEdgeColor = [0,0,1];
hs.MarkerFaceAlpha = 0.2;
hs.MarkerEdgeAlpha = 0.00;
hold on;
xx = nan(length(MCresults),1);
yy = nan(length(MCresults),1);
for t = 1:length(MCresults)
    xx(t) = MCresults{t}.ErrLowRank(k).^2;
    yy(t) = MCresults{t}.ErrLowRankEstimate{4}(k).^2;
end
hold on
hs = scatter(xx,yy);
hs.MarkerFaceColor = [1,0,0];
hs.MarkerEdgeColor = [1,0,0];
hs.MarkerFaceAlpha = 0.2;
hs.MarkerEdgeAlpha = 0.00;
plot([ErrBest(k)^2,ErrBest(k)^2],[4000,190^2]);

xlim([0^2,30000])
% set(gca,'XTick',0:25:175)

grid on;


%% Plot for k = 48
k=48;
hfig{2} = figure('Position',[100,100,350,320]);
set(hfig{2},'name','scatter-approximation-k=48','numbertitle','off');
xx = nan(length(MCresults),1);
yy = nan(length(MCresults),1);
for t = 1:length(MCresults)
    xx(t) = MCresults{t}.ErrLowRank(k).^2;
    yy(t) = MCresults{t}.ErrLowRankEstimate{2}(k).^2;
end
hold on
hs = scatter(xx,yy);
hs.MarkerFaceColor = [0,0,1];
hs.MarkerEdgeColor = [0,0,1];
hs.MarkerFaceAlpha = 0.2;
hs.MarkerEdgeAlpha = 0.00;
hold on;
xx = nan(length(MCresults),1);
yy = nan(length(MCresults),1);
for t = 1:length(MCresults)
    xx(t) = MCresults{t}.ErrLowRank(k).^2;
    yy(t) = MCresults{t}.ErrLowRankEstimate{4}(k).^2;
end
hold on
hs = scatter(xx,yy);
hs.MarkerFaceColor = [1,0,0];
hs.MarkerEdgeColor = [1,0,0];
hs.MarkerFaceAlpha = 0.2;
hs.MarkerEdgeAlpha = 0.00;
plot([ErrBest(k)^2,ErrBest(k)^2],[11^2,420]);

xlim([0^2,350])

grid on;


%% Plot for k = 128
k=128;
hfig{3} = figure('Position',[100,100,350,320]);
set(hfig{3},'name','scatter-approximation-k=128','numbertitle','off');
xx = nan(length(MCresults),1);
yy = nan(length(MCresults),1);
for t = 1:length(MCresults)
    xx(t) = MCresults{t}.ErrLowRank(k).^2;
    yy(t) = MCresults{t}.ErrLowRankEstimate{2}(k).^2;
end
hold on
hs = scatter(xx,yy);
hs.MarkerFaceColor = [0,0,1];
hs.MarkerEdgeColor = [0,0,1];
hs.MarkerFaceAlpha = 0.2;
hs.MarkerEdgeAlpha = 0.00;
hold on;
xx = nan(length(MCresults),1);
yy = nan(length(MCresults),1);
for t = 1:length(MCresults)
    xx(t) = MCresults{t}.ErrLowRank(k).^2;
    yy(t) = MCresults{t}.ErrLowRankEstimate{4}(k).^2;
end
hold on
hs = scatter(xx,yy);
hs.MarkerFaceColor = [1,0,0];
hs.MarkerEdgeColor = [1,0,0];
hs.MarkerFaceAlpha = 0.2;
hs.MarkerEdgeAlpha = 0.00;
plot([ErrBest(k)^2,ErrBest(k)^2],[0.725^2,1.125^2]);
xlim([0,1.1])
ylim([0.45,1.35])
grid on;

%%
for rr = 1:3
    set(0, 'CurrentFigure', hfig{rr})
    
    set(gca, 'FontSize', 13)
    xlabel('Approximation Error','Interpreter','latex','FontSize',17);
    ylabel('Error Estimator','Interpreter','latex','FontSize',19);
    
    set(gca,'TickDir','out')
    set(gca,'LineWidth',1,'TickLength',[0.02 0.02]);
    set(gca,'XMinorTick','on');
    set(findall(gca, 'Type', 'Line'),...
        'LineWidth',2,...
        'LineStyle','--',...
        'Color',[0.5,0.5,0.5]); % set all lines width to 2
    grid on;
    box on;
end

%% Figure TRUNCATION

%% Plot for k = 16
k = 16;
hfig{4} = figure('Position',[100,100,350,320]);
set(hfig{4},'name','scatter-truncation-k=16','numbertitle','off');
r = ceil(k/4);
% Plot for q = 5
xx = nan(length(MCresults),1);
yy = nan(length(MCresults),1);
for t = 1:length(MCresults)
    xx(t) = MCresults{t}.ErrFixedRank(k,r).^2;
    yy(t) = MCresults{t}.ErrFixedRankEstimate{2}(k,r).^2;
end
hold on
hs = scatter(xx,yy);
hs.MarkerFaceColor = [0,0,1];
hs.MarkerEdgeColor = [0,0,1];
hs.MarkerFaceAlpha = 0.2;
hs.MarkerEdgeAlpha = 0;
hs.SizeData = 48;
hold on;
% Plot for q = 10
xx = nan(length(MCresults),1);
yy = nan(length(MCresults),1);
for t = 1:length(MCresults)
    xx(t) = MCresults{t}.ErrFixedRank(k,r).^2;
    yy(t) = MCresults{t}.ErrFixedRankEstimate{4}(k,r).^2;
end
hold on
hs = scatter(xx,yy);
hs.MarkerFaceColor = [1,0,0];
hs.MarkerEdgeColor = [1,0,0];
hs.MarkerFaceAlpha = 0.2;
hs.MarkerEdgeAlpha = 0;
hs.SizeData = 48;
plot([ErrBest(r)^2,ErrBest(r)^2],[18000,72000]);

xlim([160^2,215^2])
ylim([110^2,280^2])

%% Plot for k = 48
k = 48;
hfig{5} = figure('Position',[100,100,350,320]);
set(hfig{5},'name','scatter-truncation-k=48','numbertitle','off');
r = ceil(k/4);
% Plot for q = 5
xx = nan(length(MCresults),1);
yy = nan(length(MCresults),1);
for t = 1:length(MCresults)
    xx(t) = MCresults{t}.ErrFixedRank(k,r).^2;
    yy(t) = MCresults{t}.ErrFixedRankEstimate{2}(k,r).^2;
end
hold on
hs = scatter(xx,yy);
hs.MarkerFaceColor = [0,0,1];
hs.MarkerEdgeColor = [0,0,1];
hs.MarkerFaceAlpha = 0.2;
hs.MarkerEdgeAlpha = 0;
hs.SizeData = 48;
hold on;
% Plot for q = 10
xx = nan(length(MCresults),1);
yy = nan(length(MCresults),1);
for t = 1:length(MCresults)
    xx(t) = MCresults{t}.ErrFixedRank(k,r).^2;
    yy(t) = MCresults{t}.ErrFixedRankEstimate{4}(k,r).^2;
end
hold on
hs = scatter(xx,yy);
hs.MarkerFaceColor = [1,0,0];
hs.MarkerEdgeColor = [1,0,0];
hs.MarkerFaceAlpha = 0.2;
hs.MarkerEdgeAlpha = 0;
hs.SizeData = 48;
plot([ErrBest(r)^2,ErrBest(r)^2],[900,3600]);

xlim([45.2^2,47.1^2])
ylim([500,4000])

%% Plot for k = 128
k = 128;
hfig{6} = figure('Position',[100,100,350,320]);
set(hfig{6},'name','scatter-truncation-k=128','numbertitle','off');
r = ceil(k/4);
% Plot for q = 5
xx = nan(length(MCresults),1);
yy = nan(length(MCresults),1);
for t = 1:length(MCresults)
    xx(t) = MCresults{t}.ErrFixedRank(k,r).^2;
    yy(t) = MCresults{t}.ErrFixedRankEstimate{2}(k,r).^2;
end
hold on
hs = scatter(xx,yy);
hs.MarkerFaceColor = [0,0,1];
hs.MarkerEdgeColor = [0,0,1];
hs.MarkerFaceAlpha = 0.2;
hs.MarkerEdgeAlpha = 0;
hs.SizeData = 48;
hold on;
% Plot for q = 10
xx = nan(length(MCresults),1);
yy = nan(length(MCresults),1);
for t = 1:length(MCresults)
    xx(t) = MCresults{t}.ErrFixedRank(k,r).^2;
    yy(t) = MCresults{t}.ErrFixedRankEstimate{4}(k,r).^2;
end
hold on
hs = scatter(xx,yy);
hs.MarkerFaceColor = [1,0,0];
hs.MarkerEdgeColor = [1,0,0];
hs.MarkerFaceAlpha = 0.2;
hs.MarkerEdgeAlpha = 0;
hs.SizeData = 48;
plot([ErrBest(r)^2,ErrBest(r)^2],[6^2,105]);
xlim([7.94^2,7.974^2])
ylim([5.5^2,10.5^2])

hleg = legend('$q = 5$','$q = 10$');
hleg.FontSize = 16;
hleg.Interpreter = 'latex';
hleg.Location = 'NorthWest';

hText = text(63.115,51,'Minimal Error','Interpreter','latex',...
    'FontSize',16,'Color',[0.1,0.1,0.1],'Rotation',90);

%%
for rr = 4:6
    set(0, 'CurrentFigure', hfig{rr})
    
    set(gca, 'FontSize', 13)
    xlabel('Approximation Error','Interpreter','latex','FontSize',17);
    ylabel('Error Estimator','Interpreter','latex','FontSize',19);
    
    set(gca,'TickDir','out')
    set(gca,'LineWidth',1,'TickLength',[0.02 0.02]);
    set(gca,'XMinorTick','on');
    set(findall(gca, 'Type', 'Line'),...
        'LineWidth',2,...
        'LineStyle','--',...
        'Color',[0.5,0.5,0.5]); % set all lines width to 2
    grid on;
    box on;
end


%%