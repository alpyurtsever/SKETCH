classdef Sparse < DimRedux
    %SPARSE implements a subclass of DimRedux for the sketching method
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
        function obj = Sparse(k, n, varargin)
            if ~isempty(varargin)
                if isscalar(varargin{1})
                    zeta = varargin{1};
                    varargin(1) = [];
                else
                    zeta =min(k, 8); %min(k, ceil(2*log(1 + n)));
                end
            else
                zeta = min(k, 8); % min(k, ceil(2*log(1 + n)));
            end
            if zeta<1 || zeta>n
                error('zeta should be between 1 and k.')
            end
            obj@DimRedux(k,n);
            if length(varargin) >= 1
                obj.field = lower(varargin{1});
            else
                obj.field = 'real';
            end
            
            indCol = repmat(1:n,[zeta,1]);
            indCol = indCol(:);
            indRow = nan(zeta,n);
            for tt = 1:n
                indRow(:,tt) = sort(randperm(k,zeta));
            end
            indRow = indRow(:);
            if strcmp(obj.field,'real')
                vals = sign(randn(size(indRow)));
            elseif strcmp(obj.field,'complex')
                vals = sign(randn(size(indRow)) + 1i*randn(size(indRow)));
            else
                err('Input ''field'' should be ''real'' or ''complex''.')
            end
            obj.Xi = sparse(indRow,indCol,vals,k,n,n*zeta);
            
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
            flag = true;
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

