classdef Gauss < DimRedux
    %GAUSS implements a subclass of DimRedux for the sketching method
    %described [TYUC2019]. See our paper and its supplementary materials 
    %for more details.
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
        field
        Xi
    end
    
    %% methods
    methods
        % Constructor
        function obj = Gauss(k, n, varargin)
			obj@DimRedux(k,n);
            if length(varargin) >= 1
                obj.field = lower(varargin{1});
            else
                obj.field = 'real';
            end
            if strcmp(obj.field,'real')
                obj.Xi = randn(k,n);
            elseif strcmp(obj.field,'complex')
                obj.Xi = randn(k,n) + 1i*randn(k,n);
            else
                err('Input ''field'' should be ''real'' or ''complex''.')
            end
        end
        
        %% Other methods
        function x = LeftApply(obj,A)
            x = mtimes(obj.Xi,A);
        end
            
        function x = RightApply(obj,A)
            x = mtimes(A,obj.Xi);
        end
        
        %% Overloaded methods        
        function flag = isreal(obj)
            flag = strcmp(obj.field,'real');
        end
                
        function flag = issparse(obj) %#ok
            flag = false;
        end
        
        function nz = nnz(obj)
            nz = nnz(obj.Xi);
        end
        
        function disp(obj)
            if ~istransposed(obj)
                disp(obj.Xi);
            else
                disp(obj.Xi');
            end
        end
    end
    
end

