%   This script plots the following figures in [TYUC2019]:
%   Fig.SM1
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
name = {};

%%
load data/LowRankLowNoise_n1000_R10
singVal_LowNoise = singVals;
load data/LowRankMedNoise_n1000_R10
singVal_MedNoise = singVals;
load data/LowRankHiNoise_n1000_R10
singVal_HighNoise = singVals;

name{end+1} = 'LowRankSpectra';
hfig{end+1} = figure('Position',[100,100,500,320]);
set(hfig{end},'name',name{end},'numbertitle','off');
loglog(singVal_HighNoise,'-');
hold on
loglog(singVal_MedNoise,'--');
loglog(singVal_LowNoise,'-.');

set(gca,'TickDir','out')
set(gca,'LineWidth',1,'TickLength',[0.02 0.02]);
set(findall(gca, 'Type', 'Line'),'LineWidth',2); % set all lines width to 2

ylim([.7e-5,2])

ax = gca;
ax.YTick = 10.^[-10:10];

hleg = legend('HiNoise', 'MedNoise', 'LowNoise');
hleg.FontSize = 16;
hleg.Interpreter = 'latex';
hleg.Location = 'SouthWest';

set(ax, 'FontSize', 13)
xlabel('$j$','Interpreter','latex','FontSize',17);
ylabel('$j^{th}$ singular value','Interpreter','latex','FontSize',19);

ax1 = gca;

%%
name{end+1} = 'PolyDecaySpectra';
hfig{end+1} = figure('Position',[100,100,500,320]);
set(hfig{end},'name',name{end},'numbertitle','off');
R = 10;
n = 1e3;

Sigma(1:R) = 1;
Sigma((R+1):n) = (((1:(n-R))+1).^(-0.5));
singVals_PSLOW = Sigma;
clearvars Sigma

Sigma(1:R) = 1;
Sigma((R+1):n) = (((1:(n-R))+1).^(-1));
singVals_PMED = Sigma;
clearvars Sigma

Sigma(1:R) = 1;
Sigma((R+1):n) = (((1:(n-R))+1).^(-2));
singVals_PFAST = Sigma;
clearvars Sigma

loglog(singVals_PSLOW,'-');
hold on
loglog(singVals_PMED,'--');
loglog(singVals_PFAST,'-.');

set(gca,'TickDir','out')
set(gca,'LineWidth',1,'TickLength',[0.02 0.02]);
set(findall(gca, 'Type', 'Line'),'LineWidth',2); % set all lines width to 2

hleg = legend('PolySlow', 'PolyMed', 'PolyFast');
hleg.FontSize = 16;
hleg.Interpreter = 'latex';
% hleg.Position = [0.608    0.658    0.276    0.258];
hleg.Location = 'SouthWest';
ylim([.7e-6,2])

ax = gca;
ax.YTick = 10.^[-10:1:10];

set(ax, 'FontSize', 13)
xlabel('$j$','Interpreter','latex','FontSize',17);
%ylabel('$j^{th}$ singular value','Interpreter','latex','FontSize',19);

ax.Position = ax1.Position;

%%
name{end+1} = 'ExpDecaySpectra';
hfig{end+1} = figure('Position',[100,100,500,320]);
set(hfig{end},'name',name{end},'numbertitle','off');
R = 10;
n = 1e3;

Sigma(1:R) = 1;
Sigma((R+1):n) = 10.^(-0.01*(1:(n-R)));
singVals_ESLOW = Sigma;
clearvars Sigma

Sigma(1:R) = 1;
Sigma((R+1):n) = 10.^(-0.1*(1:(n-R)));
singVals_EMED = Sigma;
clearvars Sigma

Sigma(1:R) = 1;
Sigma((R+1):n) = 10.^(-0.5*(1:(n-R)));
singVals_EFAST = Sigma;
clearvars Sigma

loglog(singVals_ESLOW,'-');
hold on
loglog(singVals_EMED,'--');
loglog(singVals_EFAST,'-.');

set(gca,'TickDir','out')
set(gca,'LineWidth',1,'TickLength',[0.02 0.02]);
set(findall(gca, 'Type', 'Line'),'LineWidth',2); % set all lines width to 2

hleg = legend('ExpSlow', 'ExpMed', 'ExpFast');
hleg.FontSize = 16;
hleg.Interpreter = 'latex';
% hleg.Position = [0.608    0.658    0.276    0.258];
hleg.Location = 'SouthWest';
ylim([.7e-9,2])
xlim([1,1e3])
ax = gca;
ax.YTick = 10.^[-10:1:10];

set(ax, 'FontSize', 13)
xlabel('$j$','Interpreter','latex','FontSize',17);
ylabel('$j^{th}$ singular value','Interpreter','latex','FontSize',19);

ax.Position = ax1.Position;


%%
name{end+1} = 'RealData';

load data/WeatherSingVals
singVals_WD = singVals;
singVals_WD = singVals_WD./max(singVals_WD);

load data/StreamVelSingVals
singVals_NS = singVals;
singVals_NS = singVals_NS./max(singVals_NS);

load data/MaxCut
singVals_MC = singVals;
singVals_MC = singVals_MC./max(singVals_MC);

load data/PhaseRetrievalSingVals
singVals_PR = diag(info.Z);
singVals_PR = singVals_PR./max(singVals_PR);


hfig{end+1} = figure('Position',[100,100,500,320]);
set(hfig{end},'name',name{end},'numbertitle','off');
loglog(singVals_NS,'-');
hold on
loglog(singVals_WD,'--');
loglog(singVals_MC,'-.');
loglog(singVals_PR,':');

set(gca,'TickDir','out');
set(gca,'LineWidth',1,'TickLength',[0.02 0.02]);
set(findall(gca, 'Type', 'Line'),'LineWidth',2); % set all lines width to 2

hleg = legend('StreamVel','MinTemp','MaxCut','PhaseRet');
hleg.FontSize = 16;
hleg.Interpreter = 'latex';
hleg.Location = 'SouthWest';

ax = gca;
ax.YTick = 10.^[-30:3:30];
xlim([1,7305])
set(ax, 'FontSize', 13)
xlabel('$j$','Interpreter','latex','FontSize',17);
%ylabel('$j^{th}$ singular value','Interpreter','latex','FontSize',19);

ax.Position = ax1.Position;

