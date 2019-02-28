classdef ThreeSketch < matlab.mixin.SetGet
    %THREESKETCH of a (m x n) dimensional target matrix A consists of three
    %matrices (X=Upsilon*A), (Y=A*Omega'), and (Z=si*A*Phi'), which are 
    %(k x n), (m x k), and (s x s) dimensional respectively. Here (k x m) 
    %dimensional matrix (Upsilon), (k x n) dimensional matrix (Omega), 
    %(s x m) dimensional matrix (Psi), and (s x n) dimensional matrix (Phi) 
    %are called as the Test Matrices. First tho sketch matrices, (X) and (Y)
    %capture information about the co-range and range of (A). The third 
    %matrix (Z) contains information about the action of (A). See our main
    %reference for more details.
    %
    %See our reference paper [TYUC2019] for the detailed explanation of the
    %sketching procedure and the arithmetic, communication and storage
    %costs.
    %
    % [TYUC2019] J.A. Tropp, A. Yurtsever, M. Udell and V. Cevher.  
    % Streaming Low-Rank Matrix Approximation with an Application to
    % Scientific Simulation. 
    %
    %Coded by: Alp Yurtsever
    %Ecole Polytechnique Federale de Lausanne, Switzerland.
    %Laboratory for Information and Inference Systems, LIONS.
    %contact: alp.yurtsever@epfl.ch
    %Created: March 01, 2018
    %Last modified: Febraury 20, 2019
    %
	% SKETCHv1.1
	% Copyright (C) 2018 Laboratory for Information and Inference Systems
	% (LIONS), Ecole Polytechnique Federale de Lausanne, Switzerland.
	% This code is a part of SKETCH toolbox. 
	% Please read COPYRIGHT before using this file.
    
    %% properties
    properties %(Access = private)
        Upsilon     % Test matrix 1
        Omega       % Test matrix 2
        Phi         % Test matrix 3
        Psi         % Test matrix 4
        Theta       % (Optional) Error test matrix
        X           % Range sketch
        Y           % Co-range sketch
        Z           % Core sketch
        W           % (Optional) Error sketch
    end
    
    %% methods
    methods
        % Constructor
        % function obj = ThreeSketch(Model, m, n, k, s, ...)
        function obj = ThreeSketch(varargin)
            Model = varargin{1};
            m = varargin{2};    % Input size: A is (m x n)
            n = varargin{3};    % Input size: A is (m x n)
            k = varargin{4};    % range parameter
            s = varargin{5};    % core parameter
            field = 'real';
            q = 0;
            for t = 6:length(varargin)
                if isnumeric(varargin{t})
                    q = varargin{t};
                elseif strcmp('real',varargin{t}) || strcmp('complex',varargin{t})
                    field = varargin{t};
                end
            end
                    
            obj.Upsilon = feval(Model,k,m,field);
            obj.Omega = feval(Model,k,n,field);
            obj.Phi = feval(Model,s,m,field);
            obj.Psi = feval(Model,s,n,field);
            
            obj.X = zeros(k,n);
            obj.Y = zeros(m,k);
            obj.Z = zeros(s,s);
            
            if q
                obj.Theta = feval('Gauss',q,m,field);
                obj.W = zeros(q,n);
            end
        end
        
        % Other methods
        LinearUpdate(obj, varargin)  
        [Q, W, P] = LowRankApprox(obj)          
        [U, S, V] = FixedRankApprox(obj, r)        
        [U, S, V] = UpaFixedRankApprox(obj, r)
        [err2] = ErrorEstimate(obj, varargin);
        [U, S, V, hf] = ScreePlot(obj, varargin);
        
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
        
        function obj = set.Z(obj,value)
            if isequal(size(value), size(obj.Z)) || isempty(obj.Z)
                obj.Z = value;
            else
                error('Size of input does not match with sketch size.')
            end
        end
        
        %% Methods about storage
        function nz = nnzTestMatrix(obj)
            nz = nnz(obj.Upsilon) + nnz(obj.Omega) ...
                + nnz(obj.Phi) + nnz(obj.Psi);
        end
        
        function nz = nnzSketch(obj)
            nz = nnz(obj.X) + nnz(obj.Y) + nnz(obj.Z);
        end
        
        function nz = nnz(obj)
            nz = nnzSketch(obj) + nnzTestMatrix(obj);
        end
        
        function nz = numelTestMatrix(obj)
            nz = numel(obj.Upsilon) + numel(obj.Omega) ...
                + numel(obj.Phi) + numel(obj.Psi);
        end
        
        function nz = numelSketch(obj)
            nz = numel(obj.X) + numel(obj.Y) + numel(obj.Z);
        end
        
        function nz = numel(obj)
            nz = numelSketch(obj) + numelTestMatrix(obj);
        end
                        
    end
end


