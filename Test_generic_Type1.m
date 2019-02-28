%Our primary goal is to obtain the best possible reconstruction of an input
%matrix using the smallest sketch possible. 
%
%For different values of storage (T), we sweep over all sketch size 
%(k,l) for sketch methods with two components, and (k,s) for sketch 
%methods with three components. Then, we record the approximation error, 
%norm(A - Aapprox,'fro'), as well as Schatten-infinity norm error for each
%case. 
%
%Test_generic_Type1 and Test_generic_Type2 files are different in the way
%input matrix A is applied. Type1 uses a dense matrix A. Type 2 computes
%the sketch of a low-rank matrix A=Hf*Hb', directly from its factors Hf
%and Hb.
%
%   See our reference paper for the detailed explanation of the
%   sketching procedure and the arithmetic, communication and storage
%   costs.
%
%	[TYUC2019] J.A. Tropp, A. Yurtsever, M. Udell and V. Cevher.  
%	Streaming Low-Rank Matrix Approximation with an Application to
%	Scientific Simulation. 
%
%   Coded by: Alp Yurtsever
%   Ecole Polytechnique Federale de Lausanne, Switzerland.
%   Laboratory for Information and Inference Systems, LIONS.
%   contact: alp.yurtsever@epfl.ch
%   Last modified: January 25, 2019
%
% SKETCHv1.1
% Copyright (C) 2018 Laboratory for Information and Inference Systems
% (LIONS), Ecole Polytechnique Federale de Lausanne, Switzerland.
% This code is a part of SKETCH toolbox. 
% Please read COPYRIGHT before using this file.


%% Beginning of the test

if strcmp(field,'real')
    alpha = 1;
elseif strcmp(field,'complex')
    alpha = 0;
else
    error('''field'' should be ''real'' or ''complex''.')
end

[m,n] = size(A);
if m < n
    A = A';
    [m,n] = size(A);
end

sd.threesketch.Tmin = (r+alpha+1)*(m+n) + (2*r+3*alpha+2)^2;
sd.twosketch.Tmin   = (r+alpha+1)*(m+n) + (alpha+1)*n;
TdivdimMin = sort([sd.threesketch.Tmin, sd.twosketch.Tmin]/(m+n));

Tdivdim = unique([r:2:2*r, round((2*r+3) * 2.^(0:.5:3.5))]);
if r == 1
    Tdivdim = unique([r:2:2*r, round((2*r+3) * 2.^(0:.5:4))]);
end
Tdivdim(Tdivdim <= TdivdimMin(1)) = [];
Tdivdim = unique([TdivdimMin,Tdivdim]);
T = Tdivdim*(m+n);

Err_S1 = (n <= 3e3);

%% Preallocation
sd.hmtsketch.l = [];
sd.hmtsketch.ErrHMT_S1 = [];
sd.hmtsketch.ErrHMT_S2 = [];
sd.hmtsketch.ErrHMT_Sinf = [];
sd.wlrtsketch.l = [];
sd.wlrtsketch.ErrWLRT_S1 = [];
sd.wlrtsketch.ErrWLRT_S2 = [];
sd.wlrtsketch.ErrWLRT_Sinf = [];

for Iter = 1:numel(Tdivdim)
        
    %% Preallocation
    sd.twosketch.k{Iter} = [];
    sd.twosketch.l{Iter} = [];
    sd.threesketch.k{Iter} = [];
    sd.threesketch.s{Iter} = [];
    sd.twosketch.ErrA7_S1{Iter} = [];
    sd.twosketch.ErrA7_S2{Iter} = [];
    sd.twosketch.ErrA7_Sinf{Iter} = [];
    sd.threesketch.ErrUpa_S1{Iter} = [];
    sd.threesketch.ErrUpa_S2{Iter} = [];
    sd.threesketch.ErrUpa_Sinf{Iter} = [];
    sd.threesketch.ErrThree_S1{Iter} = [];
    sd.threesketch.ErrThree_S2{Iter} = [];
    sd.threesketch.ErrThree_Sinf{Iter} = [];
    
    %% Sweep parameters
    k2S = r+alpha+1;
    l2S = floor((T(Iter) - k2S*m)/n);
        
    while (r+alpha+1 <= k2S) && (k2S+alpha+1 <= l2S)
        
        %% Generate sketch
        mySketch = Sketch(model, m, n, k2S, l2S, field);
        mySketch.LinearUpdate(A);
        %% Our Single-View Low-Rank Approximation (Fixed rank psd symmetric reconstruction)
        A7 = mySketch.FixedRankApprox(r);
        ErrA7_S2  = norm(A7 - A,'fro');
        if Err_S1
            Z = svd(A7 - A);
            ErrA7_S1    = norm(Z,1);
            ErrA7_Sinf  = norm(Z,inf);
        else
            ErrA7_S1    = nan;
            ErrA7_Sinf  = svds(A7 - A,1);
        end
        
        %% Print intermediate results
        fprintf('T/(m+n) = %d, k = %d, l = %d, A7RelErr = %f  ... \n', ...
            Tdivdim(Iter), k2S, l2S, ErrA7_S2./ErrBest_S2-1 );
        
        %% Save data
        sd.twosketch.k{Iter}(end+1) = k2S;
        sd.twosketch.l{Iter}(end+1) = l2S;
        sd.twosketch.ErrA7_S1{Iter}(end+1) = ErrA7_S1;
        sd.twosketch.ErrA7_S2{Iter}(end+1) = ErrA7_S2;
        sd.twosketch.ErrA7_Sinf{Iter}(end+1) = ErrA7_Sinf;
        
        %% increase k, decrease l, keep T ~constant
        k2S = k2S+1;
        l2S = floor((T(Iter) - k2S*m)/n);
        
    end
    if isempty(sd.twosketch.k{Iter})
        sd.twosketch.k{Iter} = nan;
        sd.twosketch.l{Iter} = nan;
        sd.twosketch.ErrA7_S1{Iter} = nan;
        sd.twosketch.ErrA7_S2{Iter} = nan;
        sd.twosketch.ErrA7_Sinf{Iter} = nan;
    end
    
    
    %% Three-sketch
    k3S = r+alpha+1;
    s3S = min(min(m,n), floor(sqrt(T(Iter) - k3S*(m+n))));
    if ~isreal(s3S)
        s3S = 0;
    end
    
    while (r+alpha+1 <= k3S) && (k3S+alpha <= s3S)
        
        %% Create the optimal sketch
        myThreeSketch = ThreeSketch(model, m, n, k3S, s3S, field);
        myThreeSketch.LinearUpdate(A);
        
        %% Three Single-View Fixed-Rank Approximation
        AThree = myThreeSketch.FixedRankApprox(r);
        ErrThree_S2  = norm(AThree - A,'fro');
        if Err_S1
            Z = svd(AThree - A);
            ErrThree_S1    = norm(Z,1);
            ErrThree_Sinf  = norm(Z,inf);
        else
            ErrThree_S1    = nan;
            ErrThree_Sinf  = svds(AThree - A,1);
        end
        
        %% Updahyay's Single-View Fixed-Rank Approximation
        AUpa = myThreeSketch.UpaFixedRankApprox(r);
        ErrUpa_S2  = norm(AUpa - A,'fro');
        if Err_S1
            Z = svd(AUpa - A);
            ErrUpa_S1    = norm(Z,1);
            ErrUpa_Sinf  = norm(Z,inf);
        else
            ErrUpa_S1    = nan;
            ErrUpa_Sinf  = svds(AUpa - A,1);
        end
        
        %% Compute the errors
        % Find (A - Aour) using thinSVD
        fprintf('T/(m+n) = %d, k = %d, s = %d, ThreeRelErr = %f, UpaRelErr = %f ... \n', ...
            Tdivdim(Iter), k3S, s3S, ErrThree_S2./ErrBest_S2-1, ErrUpa_S2./ErrBest_S2-1 );
        
        %% Save data
        sd.threesketch.k{Iter}(end+1) = k3S;
        sd.threesketch.s{Iter}(end+1) = s3S;
        sd.threesketch.ErrUpa_S1{Iter}(end+1) = ErrUpa_S1;
        sd.threesketch.ErrUpa_S2{Iter}(end+1) = ErrUpa_S2;
        sd.threesketch.ErrUpa_Sinf{Iter}(end+1) = ErrUpa_Sinf;
        sd.threesketch.ErrThree_S1{Iter}(end+1) = ErrThree_S1;
        sd.threesketch.ErrThree_S2{Iter}(end+1) = ErrThree_S2;
        sd.threesketch.ErrThree_Sinf{Iter}(end+1) = ErrThree_Sinf;
        
        %% increase k, decrease s, keep T ~constant
        k3S = k3S + 1;
        s3S = min(min(m,n), floor(sqrt(T(Iter) - k3S*(m+n))));
        if ~isreal(s3S)
            s3S = 0;
        end
    end
    if isempty(sd.threesketch.k{Iter})
        sd.threesketch.k{Iter} = nan;
        sd.threesketch.s{Iter} = nan;
        sd.threesketch.ErrUpa_S1{Iter} = nan;
        sd.threesketch.ErrUpa_S2{Iter} = nan;
        sd.threesketch.ErrUpa_Sinf{Iter} = nan;
        sd.threesketch.ErrThree_S1{Iter} = nan;
        sd.threesketch.ErrThree_S2{Iter} = nan;
        sd.threesketch.ErrThree_Sinf{Iter} = nan;
    end
    
    %% HMT-Sketch
    lrS = round(Tdivdim(Iter));
    sd.hmtsketch.l(end+1) = lrS;
        
    %% Generate HMT sketch
    myHmtSketch = HMTSketch(model, m, n, lrS, field);
    myHmtSketch.LinearUpdate(A);
    
    %% HMT Fixed-Rank Approximation
    AHMT = myHmtSketch.HMTFixedRankApprox(r);
    ErrHMT_S2  = norm(AHMT - A,'fro');
    if Err_S1
        Z = svd(AHMT - A);
        ErrHMT_S1    = norm(Z,1);
        ErrHMT_Sinf  = norm(Z,inf);
    else
        ErrHMT_S1    = nan;
        ErrHMT_Sinf  = svds(AHMT - A,1);
    end
    
    %% Save data
    sd.hmtsketch.ErrHMT_S1(end+1) = ErrHMT_S1;
    sd.hmtsketch.ErrHMT_S2(end+1) = ErrHMT_S2;
    sd.hmtsketch.ErrHMT_Sinf(end+1) = ErrHMT_Sinf;
    
    %% WLRT-Sketch
    sd.wlrtsketch.l(end+1) = lrS;
    
    %% Generate WLRT sketch
    myWlrtSketch = WLRTSketch(model, m, n, lrS, field);
    myWlrtSketch.LinearUpdate(A);
    
    %% WLRT Fixed-Rank Approximation
    AWLRT = myWlrtSketch.WLRTFixedRankApprox(r);
    ErrWLRT_S2  = norm(AWLRT - A,'fro');
    if Err_S1
        Z = svd(AWLRT - A);
        ErrWLRT_S1    = norm(Z,1);
        ErrWLRT_Sinf  = norm(Z,inf);
    else
        ErrWLRT_S1    = nan;
        ErrWLRT_Sinf  = svds(AWLRT - A,1);
    end
    
    %% Save data
    sd.wlrtsketch.ErrWLRT_S1(end+1) = ErrWLRT_S1;
    sd.wlrtsketch.ErrWLRT_S2(end+1) = ErrWLRT_S2;
    sd.wlrtsketch.ErrWLRT_Sinf(end+1) = ErrWLRT_Sinf;
    
    %% Print intermediate results
    fprintf('T/(m+n) = %d, l = %d, WlrtRelErr = %f, HmtRelErr = %f   ... \n', ...
        Tdivdim(Iter), lrS, ErrWLRT_S2./ErrBest_S2-1, ErrHMT_S2./ErrBest_S2-1 );
    
end

%% Variables to be saved
sd.ErrBest_S1 = ErrBest_S1;
sd.ErrBest_S2 = ErrBest_S2;
sd.ErrBest_Sinf = ErrBest_Sinf;
sd.info.m = m;
sd.info.n = n;
sd.info.r = r;
sd.info.field = field;
sd.info.model = model;
sd.info.Tdivdim = Tdivdim;
sd.info.T = T;
sd.ID = [datestr(now,30),num2str(randi([100,999]))];

