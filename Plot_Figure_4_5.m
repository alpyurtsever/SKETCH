%   This script plots the following figures in [TYUC2019]:
%   Fig.4, Fig.5
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
clearvars
close all
dname = dir('results/Test_ErrorEstimate_Figure_4_5/*.mat');
load(['results/Test_ErrorEstimate_Figure_4_5/',dname.name]);

%% Figure 4.A
hfig = {};

hfig{end+1} = figure('Position',[100,100,410,320]);
set(hfig{end},'name','tailenergy-hatA-vs-A','numbertitle','off');

spExact = sd.spectrumExact;
tau_rp1_A = sqrt(cumsum(spExact.^2,'reverse'));
h(1) = loglog(tau_rp1_A,'Color','black','LineWidth',3);
myLegend = {};
myLegend{1} = 'Actual';
hold on;

kSnap = fliplr([8,16,32,48,96,128]);
for k = kSnap
    indK = find(sd.kSweep==k);
    spApprox = sd.spectrumApprox(indK,:);
    tau_rp1_Ahat = sqrt(cumsum(spApprox.^2,'reverse'));
    h(end+1) = loglog(tau_rp1_Ahat,'LineWidth',2);       %#ok
    myLegend{end+1} = ['k = ',num2str(k)]; %#ok
end

hleg = legend(flip(h),flip(myLegend));
hleg.FontSize = 16;
hleg.Interpreter = 'latex';
hleg.Location = 'SouthWest';

ax = gca;

set(ax, 'FontSize', 13)
xlabel('Rank $(r)$','Interpreter','latex','FontSize',17);
ylabel('Tail Energy $(\tau_r)$','Interpreter','latex','FontSize',19);

set(gca,'TickDir','out')
set(gca,'LineWidth',1,'TickLength',[0.02 0.02]);

ax = gca;
ax.YTick = 10.^(-30:1:30);
ax.XTick = 10.^(-30:1:30);
xlim([1,200])
ylim([1e-3,3e3/2])
set(ax, 'FontSize', 13)
%ylabel('$j^{th}$ singular value','Interpreter','latex','FontSize',19);

grid on;
set(gca,'XMinorGrid','on');
set(gca,'YMinorGrid','off');

%% Figure 4.B

hfig{end+1} = figure('Position',[100,100,410,320]);
set(hfig{end},'name','fixed-rank-approx-error-Ar-relative','numbertitle','off');

myLegend = {};
shiftedSpectrum = circshift(sd.spectrumExact,[-1,0]);
shiftedSpectrum(end) = 0;
ErrBest = sqrt(cumsum(shiftedSpectrum.^2,'reverse'));
ErrBest = ErrBest(1:max(sd.kSweep))';
myLegend = {};
h = loglog(nan,nan);
hold on
set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
for k = kSnap
    indK = find(sd.kSweep==k);
    relErr = sd.ErrFixedRank(indK,:)./ErrBest - 1;
    loglog(relErr)
    myLegend{end+1} = ['k = ', num2str(k)];
end

% hleg = legend(myLegend);
% hleg.FontSize = 16;
% hleg.Interpreter = 'latex';
% hleg.Location = 'SouthEast';
% hleg.Position = [0.6270, 0.1531, 0.2787, 0.3803];
set(gca,'TickDir','out')
set(gca,'LineWidth',1,'TickLength',[0.02 0.02]);
set(findall(gca, 'Type', 'Line'),'LineWidth',2); % set all lines width to 2
ax = gca;
axis tight
% xlim([1,200]);

ax.YTick = 10.^(-30:1:30);
ax.XTick = 10.^(-30:1:30);
set(ax, 'FontSize', 13)
xlabel('Rank $(r)$','Interpreter','latex','FontSize',17);
ylabel('Relative Error $(S_2)$','Interpreter','latex','FontSize',19);

grid on, grid minor, grid minor;

%% Figure 5.A
%
hfig{end+1} = figure('Position',[100,100,410,320]);
set(hfig{end},'name','approx-error-hatA-error-estimate','numbertitle','off');

myColors = [0.1,0.2,1
            1,0.2,0.1
            0.8,0.4,0.2];
myLegend = {};
semilogy(sd.ErrLowRank.^2,'black','LineWidth',3);
myLegend{end+1} = 'Error in Approx.';
hold on;
% loglog(sd.ErrLowRankEstimate{2},'LineWidth',2,'Color',myColors(1,:))
% myLegend{end+1} = 'Error Est. $(q = 5)$';
loglog(sd.ErrLowRankEstimate{4}.^2,'LineWidth',2,'Color',myColors(2,:))
myLegend{end+1} = 'Error Est. $(q = 10)$';

hleg = legend(myLegend);
hleg.FontSize = 16;
hleg.Interpreter = 'latex';
hleg.Location = 'NorthEast';

set(gca,'TickDir','out')
set(gca,'LineWidth',1,'TickLength',[0.02 0.02]);
% set(findall(gca, 'Type', 'Line'),'LineWidth',2); % set all lines width to 2

ax = gca;
axis tight
ax.YTick = 10.^(-30:1:30);
% ax.XTick = 10.^(-30:1:30);
set(ax, 'FontSize', 13)
xlabel('Sketch Size $(k)$','Interpreter','latex','FontSize',17);
ylabel('Approximation Error','Interpreter','latex','FontSize',19);

grid on, grid minor, grid minor;

%% Figure 5.BCD

for k = [16,48,128]
    % First find the indices that corresponds to current choice of k
    indK = find(sd.kSweep==k);
    % Find the tail energy of A
    spExact = sd.spectrumExact;
    tau_rp1_A = sqrt(cumsum(spExact.^2,'reverse'));
    tau_rp1_A = circshift(tau_rp1_A,[-1,0]);
    tau_rp1_A(end) = 0;
    normA_fro = norm(sd.spectrumExact);
    screeActual = (tau_rp1_A./normA_fro).^2;
    screeActual = screeActual(1:k);
    % Find the tail energy of Ahat
    spApprox = sd.spectrumApprox(indK,:)';
    tau_rp1_Ahat = sqrt(cumsum(spApprox.^2,'reverse'));
    tau_rp1_Ahat = circshift(tau_rp1_Ahat,[-1,0]);
    tau_rp1_Ahat(end) = 0;
    tau_rp1_Ahat = tau_rp1_Ahat(1:k);
    % Find upper bound on tau(A)
    err2Ahat_q10 = sd.ErrLowRankEstimate{4}(indK);
    err2Ar_q10 = sd.ErrFixedRankEstimate{4}(indK,:);
    err2_0_q10 = sd.ErrZeroEstimate{4};
    screeUpper = ((tau_rp1_Ahat + err2Ahat_q10)./err2_0_q10).^2;
    % Find 
    screeLower = (tau_rp1_Ahat./err2_0_q10).^2;
    screeDirect = (err2Ar_q10./err2_0_q10).^2;
    % Create the figure handles
    hfig{end+1} = figure('Position',[100,100,410,320]);
    set(hfig{end},'name',['scree-plots-k=',num2str(k)],'numbertitle','off');
    % Plot lines
    h1 = semilogy(screeActual,'Color','black','LineWidth',3);
    hold on
    h2 = semilogy(screeLower,'--','Color',[1,0.1,0.1],'LineWidth',2);
    h3 = semilogy(screeUpper,':','Color','blue','LineWidth',2);
%     semilogy(screeDirect,'Color',[0.7,0.6,0],'LineWidth',2);
    if k == 128
    hleg = legend([h3,h2,h1],...
         'Upper (6.10)',...
         'Lower (6.9)',...
    	 'Actual (6.8)'...
         );
    hleg.FontSize = 16;
    hleg.Interpreter = 'latex';
    hleg.Location = 'SouthWest';
    end
    % Set axis properties
    ax = gca;
    set(ax,'TickDir','out')
    set(ax,'LineWidth',1,'TickLength',[0.02 0.02]);
    axis tight
    ax.YTick = 10.^(-30:1:30);
    set(ax, 'FontSize', 13)
    xlabel('Rank $(r)$','Interpreter','latex','FontSize',17);
    ylabel('Proportion of Energy','Interpreter','latex','FontSize',19);
    grid on, grid minor, grid minor;
end
