classdef RVEData < handle
    properties
        name                 % Name of the RVE
        phases               % Name of phases homogenized
        homorder             % Order of homogenization
        shape                % shape and orientation
        tol                  % precision of computation
        C0
        %P                    % Eshelby problem-related strain concentration tensors
        %Ainf                 % Eshelby problem-related strain concentration tensors
        %EEinfty              % link tensor
        %A                    % phase strain concentration tensors
        %APstroud             % if stroud integration necessary, product of A and P
        %C                    % homogenized stiffness tensor
    end
    properties (Dependent)
        vol                  % volume fractions of RVE phase
        f                    % phase volume fractions related to RVE       
    end
    
    methods
        %% CONSTRUCTOR
        function newRVE = RVEData(namehom,phaseshom,homorderhom,shapehom,toli)
            if nargin > 0
                newRVE.name = namehom;
                % validation for phase input
                newRVE.phases=[];
                for i=1:numel(phaseshom)
                    if isa(phaseshom(i),'PhaseData')
                        newRVE.phases.(phaseshom(i).name) = phaseshom(i);
                    else
                        error('Phase must be a previously defined Material object, use PhaseData for definition')
                    end
                end
                newRVE.homorder = homorderhom;
                % validation input for shape
                valid_shape={'sphere','needle','disc','spheroid','ellipsoid','matrix'};
                if ismember(shapehom{1},valid_shape)
                    newRVE.shape = shapehom;
                end
                % validation for tolerance input 
                if toli>=1 && toli<=16
                newRVE.tol=toli;
                else
                    error('Tolerance must be between 1 (very low) and 16 (very high)')
                end
            end
        end
        %% VOLUME FRACTIONS
        function vol = get.vol(obj)
            vol=0;
            fld=fieldnames(obj.phases);
            % sum up all phase volumes
            for i=1:numel(fld)
                vol=vol+obj.phases.(fld{i}).vol;
            end
        end

        function f = get.f(obj)
            fld=fieldnames(obj.phases);
            % get RVE-related phase volume fractions
            for i=1:numel(fld)
                f.(fld{i})=obj.phases.(fld{i}).vol/obj.vol;
            end
        end
        function obj = set.f(obj,~)
            fprintf('%s%d\n','volumes read as: ',obj.f)
            error('You cannot set volumes of RVE, use phase volume definitions');
        end
        
        
        %% HOMOGENIZATION
        function [phasehom,C,A,Ainf,P,APstroud,EEinfty] = homE(obj,scheme)
            %gengerate string and function handle for homogenization scheme
            str=['fun_hom_',scheme];
            fhom=str2func(str);           
            %evaluate homogenization
            [C,A,Ainf,P,APstroud,EEinfty,C0]=fhom(obj);      
            
            % generate new material and phase
            if sum(sum(abs(C-transpose(C))))>1e-10
                error('homogenized stiffness tensor is not symmetric - something went wrong')
            end
            Mathom=MatData(obj.name,'hom',(C+transpose(C))/2,false,false);
            phasehom=PhaseData(obj.name,Mathom,obj.shape,obj.vol);
        end
        function [phasehom,C,A,Ainf,P,APstroud,EEinfty] = eigstraininfl(obj,scheme)           
            %gengerate string and function handle for homogenization scheme
            str=['fun_eigstraininfl_',scheme];
            fhom=str2func(str);
            %evaluate homogenization
            [C,A,Ainf,P,APstroud,EEinfty]=fhom(obj);      
            
            % generate new material and phase
            if sum(sum(abs(C-transpose(C))))>1e-10
                error('homogenized stiffness tensor is not symmetric - something went wrong')
            end
            Mathom=MatData(obj.name,'hom',(C+transpose(C))/2,false,false);
            phasehom=PhaseData(obj.name,Mathom,obj.shape,obj.vol);
        end
        
    end
  
end