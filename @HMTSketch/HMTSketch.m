classdef HMTSketch < matlab.mixin.SetGet
    %HMTSKETCH of a (m x n) dimensional target matrix A takes the form
    %(Y = A*Omega) and (Ytilde = A'*OmegaTilde) which are (n x l) and (m x l)
    %dimensional respectively. Here (m x l) dimensional matrix (Omega) and
    %(n x l) dimensional matrix (OmegaTilde) are called as the Test Matrices.
    %This code implements methods from [HMT2011]. 
    %
    %[HMT2011] N. Halko, P.G. Martinsson and J.A. Tropp. Finding Structure
    %with Randomness: Probabilistic algorithms for constructing approximate 
    %matrix decompositions,
    %
    %[TYUC2019] J.A. Tropp, A. Yurtsever, M. Udell and V. Cevher.  
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
    %SKETCHv1.1
	%Copyright (C) 2018 Laboratory for Information and Inference Systems
	%(LIONS), Ecole Polytechnique Federale de Lausanne, Switzerland.
	%This code is a part of SKETCH toolbox. 
	%Please read COPYRIGHT before using this file.

    %% properties
    properties (Access = private)
        Omega           % (m x l) dimensional test matrix
        OmegaTilde      % (n x l) dimensional test matrix
        Y           % (n x l) dimensional sketch
        Ytilde      % (m x l) dimensional sketch
    end
    
    %% methods
    methods
        % Constructor
        function obj = HMTSketch(varargin)
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
            obj.Omega = feval(Model,l,n,field)';
            obj.OmegaTilde = feval(Model,l,m,field)';
            obj.Y = zeros(m,l);
            obj.Ytilde = zeros(n,l);
        end
        
        % Other methods
        LinearUpdate(obj, varargin)                 % Algorithm 2 in [MR]
        
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
            nz = nnz(obj.Omega) + nnz(obj.Omegatilde);
        end
        
        function nz = nnzSketch(obj)
            nz = nnz(obj.Y) + nnz(obj.Ytilde);
        end
        
        function nz = nnz(obj)
            nz = nnzSketch(obj) + nnzTestMatrix(obj);
        end
        
        function nz = numelTestMatrix(obj)
            nz = numel(obj.Omega) + numel(obj.Omegatilde);
        end
        
        function nz = numelSketch(obj)
            nz = numel(obj.Y) + numel(obj.Ytilde);
        end
        
        function nz = numel(obj)
            nz = numelSketch(obj) + numelTestMatrix(obj);
        end
    end
end


