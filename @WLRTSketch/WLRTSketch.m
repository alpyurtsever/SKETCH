classdef WLRTSketch < matlab.mixin.SetGet
    %WLRTSKEKTCH of a (m x n) dimensional target matrix A takes the form
    %(Y = R*A) and (Ytilde = Rtilde*A') which are (l x n) and (l x m)
    %dimensional respectively. Here (l x m) dimensional matrix (R) and
    %(l x n) dimensional matrix (Rtilde) are called as the Test Matrices.
    %This code implements methods from [WLRT2008]. 
    %
    %[WLRT2008] F. Woolfe, E. Liberty, V. Rokhlin, M. Tygert. A Fast
    %Randomized Algorithm for the Approximation of Matrices.
    %
	%[TYUC2019] J.A. Tropp, A. Yurtsever, M. Udell and V. Cevher, 
    %Streaming Low-Rank Matrix Approximation with an Application to 
    %Scientific Simulation.
    %
    %Coded by: Alp Yurtsever
    %Ecole Polytechnique Federale de Lausanne, Switzerland.
    %Laboratory for Information and Inference Systems, LIONS.
    %contact: alp.yurtsever@epfl.ch
    %Created: January 18, 2019
    %Last modified: January 18, 2019
    %
	% SKETCHv1.1
	% Copyright (C) 2018 Laboratory for Information and Inference Systems
	% (LIONS), Ecole Polytechnique Federale de Lausanne, Switzerland.
	% This code is a part of SKETCH toolbox. 
	% Please read COPYRIGHT before using this file.
    
    %% properties
    properties (Access = private)
        R           % (l x m) dimensional test matrix (R = Omega')
        Rtilde      % (l x n) dimensional test matrix (Rtilde = OmegaTilde')
        Y           % (l x n) dimensional sketch
        Ytilde      % (l x m) dimensional sketch
    end
    
    %% methods
    methods
        % Constructor
        function obj = WLRTSketch(varargin)
            narginchk(4,5);
            Model = varargin{1};
            m = varargin{2};
            n = varargin{3};
            l = varargin{4};
            if nargin == 5
                field = varargin{5};
            else
                field = 'real';
            end
            obj.R = feval(Model,l,m,field);
            obj.Rtilde = feval(Model,l,n,field);
            obj.Y = zeros(l,n);
            obj.Ytilde = zeros(l,m);
        end
        
        % Other methods
        LinearUpdate(obj, varargin)
        [U, Sigma, V] = WLRTFixedRankApprox(obj, r)
        
        %% Property set methods
        function obj = set.Ytilde(obj,value)
            if isequal(size(value), size(obj.Ytilde)) || isempty(obj.Ytilde)
                obj.Ytilde = value;
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
            nz = nnz(obj.R) + nnz(obj.Rtilde);
        end
        
        function nz = nnzSketch(obj)
            nz = nnz(obj.Y) + nnz(obj.Ytilde);
        end
        
        function nz = nnz(obj)
            nz = nnzSketch(obj) + nnzTestMatrix(obj);
        end
        
        function nz = numelTestMatrix(obj)
            nz = numel(obj.R) + numel(obj.Rtilde);
        end
        
        function nz = numelSketch(obj)
            nz = numel(obj.Y) + numel(obj.Ytilde);
        end
        
        function nz = numel(obj)
            nz = numelSketch(obj) + numelTestMatrix(obj);
        end
    end
end


