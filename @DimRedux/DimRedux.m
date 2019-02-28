classdef DimRedux
    %DIMREDUX implements a class definition for the dimensionality  
    %reduction maps described in [TYUC2019]. See our paper and its
    %supplementary materials for more details. 
    %
    %[TYUC2019] J.A. Tropp, A. Yurtsever, M. Udell and V. Cevher.  
    %Streaming Low-Rank Matrix Approximation with an Application to
    %Scientific Simulation. 
    %
    %Coded by: Alp Yurtsever
    %Ecole Polytechnique Federale de Lausanne, Switzerland.
    %Laboratory for Information and Inference Systems, LIONS.
    %contact: alp.yurtsever@epfl.ch
    %Last modified: July 22, 2018
    %
    %SKETCHv1.1
	%Copyright (C) 2018 Laboratory for Information and Inference Systems
	%(LIONS), Ecole Polytechnique Federale de Lausanne, Switzerland.
	%This code is a part of SKETCH toolbox. 
	%Please read COPYRIGHT before using this file.
    
    %% properties
    properties (Access = private)
        k
        n
        transposeFlag
    end
    
    %% methods
    methods
        % Constructor
        function obj = DimRedux(k, n) 
            if k > n
                error('k should be less than or equal to n.')
            end
            obj.k = k;
            obj.n = n;
            obj.transposeFlag = 0;
        end
        
        %% Mtimes
        function x = mtimes(obj1,obj2)
            if isa(obj1,'DimRedux') && isa(obj2,'double')
                if ~obj1.transposeFlag
                    x = LeftApply(obj1,obj2);
                else
                    x = RightApply(obj1',obj2')';
                end
            elseif isa(obj1,'double') && isa(obj2,'DimRedux')
                if ~obj2.transposeFlag
                    x = RightApply(obj2,obj1);
                else
                    x = LeftApply(obj2',obj1')';
                end
            else
                error(['Multiplication is defined only between', ...
                    'DimRedux and double classes.']);
            end
        end
        
        %% Plus 
        % (this can be implemented computation and storage efficiently)
        function x = plus(obj1,obj2)
            if isa(obj1,'DimRedux') && isa(obj2,'double')
                if obj1.transposeFlag
                    A = obj1*speye(size(obj1,2));
                else
                    A = speye(size(obj1,1))*obj1;
                end
                x = A + obj2;
            elseif isa(obj1,'double') && isa(obj2,'DimRedux')
                if obj2.transposeFlag
                    A = obj2*speye(size(obj2,2));
                else
                    A = speye(size(obj2,1))*obj2;
                end
                x = A + obj1;
            else
                error(['Addition is defined only between', ...
                    'DimRedux and double classes.']);
            end
        end
        
        %% return a flag for transposition
        function flag = istransposed(obj)
            flag = obj.transposeFlag;
        end
                
        %% Overloaded methods
        function obj = ctranspose(obj)
            obj.transposeFlag = ~obj.transposeFlag;
        end
        
        function [k,n] = size(obj,varargin)
            if ~obj.transposeFlag
                k = obj.k;
                n = obj.n;
            else
                k = obj.n;
                n = obj.k;
            end
            if nargin == 2
                nargoutchk(0,1)
                if ~isempty(varargin)
                    if varargin{1} == 1
                        n = [];
                    elseif varargin{1} == 2
                        k = n;
                        n = [];
                    else
                        error('Second argument should be 1 or 2.');
                    end
                end
            elseif nargout <= 1
                k = [k,n];
            end
        end
        
        function nm = numel(obj)
            nm = obj.n*obj.k;
        end
        
        function ln = length(obj)
            ln = obj.n; % k cannot be larger than n ==> max(size(obj)) = n
        end
    end
end

