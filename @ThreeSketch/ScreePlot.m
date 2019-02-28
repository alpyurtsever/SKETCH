function [U,Sigma,V, hf] = ScreePlot( obj, varargin )
%SCREEPLOT plots the upper and lower estimates of the scree plot. Then,
%upon the rank input from the user, it returns the fixed rank approximation
%from the sketch in the factorized form, U*Sigma*V'. See Section 6.5 from 
%our main reference [TYUC2019] for more details. 
%
%   Require: (optional) truncation rank r
%   Ensure:  Returns factors 
%             (m x r) dimensional U with orthonormal columns
%             (r x r) dimensional nonnegative diagonal matrix Sigma
%             (n x r) dimensional V with orthonormal columns
%            that form a rank-r approximation (Aout = U*Sigma*V') of the 
%            sketched matrix.
%
%   [U,Sigma,V] = S.ScreePlot() draws the scree plot for sketch S, and
%   waits for the users keyboard input for the truncation rank. Then, it
%   returns the factors U, Sigma, V that forms (Aout = U*Sigma*V').
%
%   [U,Sigma,V,hf] = S.ScreePlot() is the same routine, also returns the
%   handle hf for scree plot. 
%
%   Aout = S.ScreePlot() is the same routine, but returns Aout in the
%   ambient dimensions, instead of the factors. [Aout,hf]=S.ScreePlot()
%   also returns the figure handle. 
%
%   For all cases above, S.ScreePlot(r) sets the truncation rank in
%   advance, hence surpases for the user keyboard input. 
%   
% [TYUC2019] J.A. Tropp, A. Yurtsever, M. Udell and V. Cevher.  
% Streaming Low-Rank Matrix Approximation with an Application to
% Scientific Simulation. 
%
%Coded by: Alp Yurtsever
%Ecole Polytechnique Federale de Lausanne, Switzerland.
%Laboratory for Information and Inference Systems, LIONS.
%contact: alp.yurtsever@epfl.ch
%Created: February 05, 2019
%Last modified: February 22, 2019
%
% SKETCHv1.1
% Copyright (C) 2018 Laboratory for Information and Inference Systems
% (LIONS), Ecole Polytechnique Federale de Lausanne, Switzerland.
% This code is a part of SKETCH toolbox. 
% Please read COPYRIGHT before using this file.

narginchk(1,2);

if isreal(obj.Theta)
    beta = 1;
else
    beta = 2;
end

% Get Sketch parameters
q = size(obj.W,1);
if ~q
    error('You need to draw Error Sketch to use ScreePlot feature.');
end

% Compute low-rank approximation
[Q,W,P] = obj.LowRankApprox();

% Compute error estimate
err2 = obj.ErrorEstimate(Q,W,P);
err2_0 = obj.ErrorEstimate(0);

% Use rank truncation
[U, Sigma, V] = svd(W);

% Compute error estimate after truncation
err2_r = nan(size(obj.X,1),1);
for r = 1:size(obj.X,1)
    UU = Q*U(:,1:r);
    SS = Sigma(1:r,1:r);
    VV = P*V(:,1:r);
    err2_r(r) = obj.ErrorEstimate(UU,SS,VV);
end

% Compute tail energy
tau_rp1_Ahat = tailEnergy(diag(Sigma));

% Compute scree error (Uppet)
screeUpper = ((tau_rp1_Ahat + err2)/err2_0).^2;

% Compute scree error (Lower)
screeLower = (tau_rp1_Ahat/err2_0).^2;

% Compute scree error (Approx)
% screeDirect = (err2_r/err2_0).^2;

% Draw ScreePlot
hf = figure('Position',[100,100,550,220]);
set(hf,'name','ScreePlot-SST','numbertitle','off');
h1 = plot(screeUpper,':','LineWidth',2,'Color','blue');
hold on;
h2 = plot(screeLower,'--','LineWidth',2,'Color',[1,0.2,0.1]);
% h3 = semilogy(screeDirect,'-','LineWidth',1.5,'Color',[0.3,0.3,0.3]);
ax = gca;
set(ax,'TickDir','out')
% set(ax,'YTick',10.^(-100:100));
set(ax,'LineWidth',1,'TickLength',[0.02 0.02]);
set(ax, 'FontSize', 13)
xlabel('Rank $(r)$','Interpreter','latex','FontSize',17);
ylabel('Proportion of Energy','Interpreter','latex','FontSize',16);
grid on, grid minor, grid minor;
axis tight;
hleg = legend([h2,h1], 'Lower (6.9)', 'Upper (6.10)');
hleg.FontSize = 16;
hleg.Interpreter = 'latex';
hleg.Location = 'NorthEast';
% Choose rank r
if isempty(varargin)
    % Get the rank-input
    set(hf,'visible','on');
    drawnow;
    r = input('Input the rank (find a knee in ScreePlot): ');
else
    r = varargin{1};
end
hline = plot([1,1]*r, ylim,'LineWidth',1,'Color',[0.4,0.4,0.4]);
uistack(hline,'bottom');

% Truncate and output
Sigma = Sigma(1:r,1:r);
U = Q*U(:,1:r);
V = P*V(:,1:r);

% Output
if nargout == 1; U = U*Sigma*V'; end
if nargout == 2; U = U*Sigma*V'; Sigma = hf; end

end

function tau = tailEnergy(spectrum)
    tau = sqrt(cumsum(spectrum.^2,'reverse'));
    tau(1:end-1) = tau(2:end); % we can use circshift (faster but less readable)
    tau(end) = 0;
end