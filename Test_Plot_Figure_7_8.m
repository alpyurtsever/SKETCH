%   This script implements the experiments, and plots the following figures 
%   in [TYUC2019]:
%   Fig.7, Fig.8
%
%   DATA LINK: Please download the "daily mean" and "land mask" data from 
%   https://www.esrl.noaa.gov/psd/data/gridded/data.noaa.oisst.v2.highres.html
%   into the data folder.
%
%	[TYUC2019] J.A. Tropp, A. Yurtsever, M. Udell and V. Cevher.  
%	Streaming Low-Rank Matrix Approximation with an Application to
%	Scientific Simulation. 
%
%   Coded by: Alp Yurtsever
%   Ecole Polytechnique Federale de Lausanne, Switzerland.
%   Laboratory for Information and Inference Systems, LIONS.
%   contact: alp.yurtsever@epfl.ch
%   Created: February 5, 2019
%   Last modified: February 22, 2019
%
% SKETCHv1.1
% Copyright (C) 2018 Laboratory for Information and Inference Systems
% (LIONS), Ecole Polytechnique Federale de Lausanne, Switzerland.
% This code is a part of SKETCH toolbox. 
% Please read COPYRIGHT before using this file.

% Preamble
rng(0,'twister');
clearvars;
close all;
load utils/icyhot;  % Colormap for geoplots

%% Load or Create and Save data
if exist('data/SST-sketch.mat','file')
    load('data/SST-sketch.mat');
else
    % Consider years from 1981 to 2019
    years = 1981:2019;
    
    % Data dimensions (m x n)
    m = 691150; % points on the map
    n = 13670;  % time instance (days)
    
    % Choose sketch-size parameters
    alpha = 1; % since we are using real field
    Tdivdim = 48;
    T = Tdivdim*(m+n);
    
    k3S = floor( (1/8) * (sqrt((m+n+4*alpha)^2 + 16*(T - alpha^2)) - (m+n+4*alpha)) );
    s3S = floor(sqrt(T - k3S*(m+n)));
    q3S = 10;
    
    % Draw test matrices
    myThreeSketch = ThreeSketch('Sparse', m, n, k3S, s3S, q3S, 'real');
    
    % Load the mask to omit lands
    lsmask = find(ncread('data/lsmask.oisst.v2.nc','lsmask'));
    
    % Load data (stream as years)
    day = 0; % Total number of days so far
    day2year(1) = 0;
    for yr = years
        
        fprintf('Streaming year %d. \n',yr);
        filename = ['data/sst.day.mean.',num2str(yr),'.nc'];
        % time = ncread(filename,'time');
        sst = ncread(filename,'sst');
        dataStream = sparse([],[],[],m,n,m*size(sst,3));
        for t = 1:size(sst,3)
            day = day + 1;
            sst_dy = sst(:,:,t);
            dataStream(:,day) = sst_dy(lsmask);  %#ok: can be done much more efficiently but doesn't worth the trouble!
        end
        day2year(end+1) = day; %#ok
        myThreeSketch.LinearUpdate(dataStream,1,1);
        clearvars dataStream;
        clearvars sst_dy;
        clearvars sst;
        
    end
    
    % Load latitude and longitude (we will need to plot afterwards)
    lat = ncread(filename,'lat');
    lon = ncread(filename,'lon');
    clearvars filename;
    clearvars yr dy;
    
    % Save the sketch
    save('data/SST-sketch.mat','-v7.3');
end


%% Use ScreePlot to approximate rank
[U,Sigma,V,hf] = myThreeSketch.ScreePlot(5);

%% Create Figures with Geoshow Plotting

hfig = {};

latitude = repmat(lat',[1440,1]);
longitude = repmat(lon,[1,720]);

%% Plot left singulat vectors (geoplots)
for r = 1:size(U,2)
    Z = nan(1440,720);
    Z(lsmask) = -U(:,r);
    
    hfig{r} = figure; %#ok
    set(hfig{r},'name',['SST-spatial-r=',num2str(r)],'numbertitle','off');
    axesm eckert4;
    framem; gridm;
    axis off
    
    geoshow(latitude,longitude,Z,'DisplayType','texturemap');
    geoshow('landareas.shp', 'FaceColor', 'white');
    colormap(icyhot);
    cbarlimit = max(abs(Z(:)));
    caxis([-1, 1]*cbarlimit);
end

%% Plot right singular vectors (time series)
for r = 1:size(V,2)
    
    hfig{r} = figure('Position',[100,100,600,250]); %#ok
    set(hfig{r},'name',['SST-time-r=',num2str(r)],'numbertitle','off');
    
    if r > 1
    hline = plot([1,n], [0,0],'LineWidth',1.5,'Color',[0.4,0.4,0.4]);
    hold on
    end
    plot(-V(:,r),'LineWidth',2,'Color','black');

    axis tight
    grid on; grid minor; grid minor
    
    ax = gca;
    
    yearsSample = 1985:5:2019;
    indSample = [];
    for tt = 1:length(yearsSample)
        indSample(end+1) = find(years == yearsSample(tt)); %#ok
    end
    daysSample = day2year(indSample);
    set(ax,'XTick',daysSample,'XTickLabel',yearsSample);
    set(ax,'TickDir','out')
    set(ax,'LineWidth',1,'TickLength',[0.02 0.02]);    
    set(ax, 'FontSize', 13)
    set(ax,'TickDir','out')
    set(ax,'LineWidth',1,'TickLength',[0.02 0.02]);
    set(ax,'XMinorTick','on');
    grid on;
    box on;
    
    if r == 1
        ylims = ylim;
        ylims = linspace(ylims(1), ylims(2) ,5);
        set(ax,'YTick',ylims)
    else
        ylims = max(abs(ylim))*[-1,1];
        ylim(ylims);
        ylims = linspace(ylims(1), ylims(2) ,5);
        set(ax,'YTick',ylims)
    end
    
end

%% Table
axObjs = hf.Children; %Scree Plot
lines = findall(axObjs, 'Type', 'Line'); % Get Lines From Scree Plot
screeLower = lines(1).YData';
screeUpper = lines(2).YData';
screeLower = screeLower(1:12);
screeUpper = screeUpper(1:12);
[screeLower, screeUpper]

