classdef NystromSketch < matlab.mixin.SetGet
    %NYSTROMSKETCH implements a class definition for the sketching method 
    %described [TYUC2017Nys]. 
    %
	%[TYUC2017Nys] J.A. Tropp, A. Yurtsever, M. Udell and V. Cevher. Fixed-
	%Rank Approximation of a Positive-Semidefinite Matrix from Streaming 
	%Data. In Proc. 31st Conference on Neural Information Processing Systems
	%(NeurIPS), Long Beach, CA, USA, December 2017.
    %
    %Coded by: Alp Yurtsever
    %Ecole Polytechnique Federale de Lausanne, Switzerland.
    %Laboratory for Information and Inference Systems, LIONS.
    %contact: alp.yurtsever@epfl.ch
    %Created: April 12, 2017
    %Last modified: July 16, 2018 (Modified for SKETCHv1.0, by Alp Yurtsever)
    %
    %NysSKETCHv1.0
	%Copyright (C) 2017 Laboratory for Information and Inference Systems
	%(LIONS), Ecole Polytechnique Federale de Lausanne, Switzerland.
	%This code is a part of Nys???SKETCH toolbox. 
	%Please read COPYRIGHT before using this file.
    %This code is a part of SKETCH toolbox. 
    %Please read COPYRIGHT before using this file.
	

    
    %% properties
    properties (Access = private)
        Omega     % (n x k) dimensional test matrix for the range of A (std Gaussian + orthonormalization)
        Y         % (n x k) dimensional range sketch
    end
    
    %% methods
    methods
        % Constructor
        % function obj = NystromSketch(A, k, varargin)
        % function obj = Sketch(Model, n, k, ...)
        function obj = NystromSketch(varargin)
            narginchk(3,4);
            Model = varargin{1};
            n = varargin{2};
            k = varargin{3};
            t = 4;
            if nargin == (t)
                field = varargin{t};
            else
                field = 'real';
            end
            obj.Omega = feval(Model,k,n,field);
            obj.Omega = obj.Omega';
            obj.Y = zeros(n,k);
        end
        
        % Other methods
        LinearUpdate(obj, varargin)
        [Q, D]  = FixedRankPSDApprox(obj, r)
        [Q, D]  = GittensMahoneyApprox(obj, r)
        [Q, D]  = FixedRankPSDApproxUnstable(obj, r)
        
        %% Property set methods
        function obj = set.Y(obj,value)
            if isequal(size(value), size(obj.Y)) || isempty(obj.Y)
                obj.Y = value;
            else
                error('Size of input does not match with sketch size.')
            end
        end
        
        %% Methods about storage
        function nz = nnzTestMatrix(obj)
            nz = nnz(obj.Omega);
        end
        
        function nz = nnzSketch(obj)
            nz = nnz(obj.Y);
        end
        
        function nz = nnz(obj)
            nz = nnzSketch(obj) + nnzTestMatrix(obj);
        end
        
        function nz = numelTestMatrix(obj)
            nz = numel(obj.Omega);
        end
        
        function nz = numelSketch(obj)
            nz = numel(obj.Y);
        end
        
        function nz = numel(obj)
            nz = numelSketch(obj) + numelTestMatrix(obj);
        end
        
    end
end

