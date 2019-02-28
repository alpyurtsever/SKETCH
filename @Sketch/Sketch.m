classdef Sketch < matlab.mixin.SetGet
    %SKETCH of a (m x n) dimensional target matrix A takes the form
    %(Y = A*Omega) and (X = Upsilon*A) which are (m x k) and (l x n)
    %dimensional respectively. Here (n x k) dimensional matrix (Omega) and
    %(l x m) dimensional matrix (Upsilon) are called as the Test Matrices. This
    %codes implements the sketching method proposed in the reference paper,
    %where the test matrices are statistically independent and follow the
    %standard normal distribution. The matrix (Y) collects information
    %about the range of (A), while the matrix (X) collects information
    %about the co-range of (A). Both parts are necessary.
    %
	%   NOTICE: THIS FILE IS MODIFIED FROM THE PRACTICALSKETCHING TOOLBOX 
	%   by Alp Yurtsever to add SSFT test matrices. Nameing conventions are
	%   also changed to match with the new paper.
    %
    %See our reference paper [TYUC2017P] for the detailed explanation of the
    %sketching procedure and the arithmetic, communication and storage
    %costs.
    %
	%[TYUC2017P] J.A. Tropp, A. Yurtsever, M. Udell and V. Cevher, Practical  
	%Sketching Algorithms for Low-Rank Matrix Approximation, SIMAX, 2017.
    %
    %Coded by: Alp Yurtsever
    %Ecole Polytechnique Federale de Lausanne, Switzerland.
    %Laboratory for Information and Inference Systems, LIONS.
    %contact: alp.yurtsever@epfl.ch
    %Created: June 28, 2016
    %Last modified: July 16, 2018 (Modified for SKETCHv1.0, by Alp Yurtsever)
    %
    %PRACTICALSKETCHING-v1.0
    %Copyright (C) 2017 Laboratory for Information and Inference Systems
    %(LIONS), Ecole Polytechnique Federale de Lausanne, Switzerland.
    %
    %This code is a part of PRACTICALSKETCHING toolbox.
    %Please read COPYRIGHT before using this file.
    %This code is a part of SKETCH toolbox. 
    %Please read COPYRIGHT before using this file.

    %% properties
    properties (Access = private)
        Omega     % (k x n) dimensional test matrix for the range of A
        Upsilon   % (l x m) dimensional test matrix for the corange of A
        Y         % (m x k) dimensional range sketch
        X         % (l x n) dimensional corange sketch
    end
    
    %% methods
    methods
        % Constructor
        % function obj = Sketch(A, k, l, Field, Orthogonalization)
%         function obj = Sketch(A, k, l, varargin)
        % function obj = Sketch(Model, m, n, k, l, ...)
        function obj = Sketch(varargin)
            narginchk(5,6);
            Model = varargin{1};
            m = varargin{2};
            n = varargin{3};
            k = varargin{4};
            l = varargin{5};
            if nargin == 6
                field = varargin{6};
            else
                field = 'real';
            end
            obj.Omega = feval(Model,k,n,field);
            obj.Upsilon = feval(Model,l,m,field);
            obj.Y = zeros(m,k);
            obj.X = zeros(l,n);
        end
        
        % Other methods
        LinearUpdate(obj, varargin)                 % Algorithm 2 in [MR]
        [Q, W] = SimpleLowRankApprox(obj)           % Algorithm 3
        [Q, W] = LowRankApprox(obj)                 % Algorithm 4
        [U, S] = LowRankSymApprox(obj)              % Algorithm 5
        [U, D] = LowRankPSDApprox(obj)              % Algorithm 6
        [Q, Sigma, V] = FixedRankApprox(obj, r)     % Algorithm 7
        [U, D] = FixedRankSymApprox(obj, r)         % Algorithm 8
        [U, D] = FixedRankPSDApprox(obj, r)         % Algorithm 9
        
        %% Property set methods
        function obj = set.X(obj,value)
            if isequal(size(value), size(obj.X)) || isempty(obj.X)
                obj.X = value;
            else
                error('Size of input does not match with sketch size.')
            end
        end
        
        function obj = set.Y(obj,value)
            if isequal(size(value), size(obj.Y)) || isempty(obj.Y)
                obj.Y = value;
            else
                error('Size of input does not match with sketch size.')
            end
        end
        
        %% Methods about storage
        function nz = nnzTestMatrix(obj)
            nz = nnz(obj.Omega) + nnz(obj.Upsilon);
        end
        
        function nz = nnzSketch(obj)
            nz = nnz(obj.Y) + nnz(obj.X);
        end
        
        function nz = nnz(obj)
            nz = nnzSketch(obj) + nnzTestMatrix(obj);
        end
        
        function nz = numelTestMatrix(obj)
            nz = numel(obj.Omega) + numel(obj.Upsilon);
        end
        
        function nz = numelSketch(obj)
            nz = numel(obj.Y) + numel(obj.X);
        end
        
        function nz = numel(obj)
            nz = numelSketch(obj) + numelTestMatrix(obj);
        end
    end
end


