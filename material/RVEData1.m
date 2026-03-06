classdef RVEData1 < handle
    properties
        name                 % Name of the RVE
        phases               % Name of phases homogenized
        scheme               % Homogenization scheme
        %homorder            % Order of homogenization
        shape                % shape and orientation of the homogenized phase
        tol                  % precision of computation
        P                    % Eshelby problem-related strain concentration tensors
        Ainf                 % Eshelby problem-related strain concentration tensors
        EEinfty              % link tensor
        A                    % phase strain concentration tensors
        APstroud             % if stroud integration necessary, product of A and P
        C0                   % stiffness of infinite matrix in matrix-inclusion problem
        Chom                 % homogenized stiffness tensor
        Q                    % eigenstrain influence tensors of phases
        PIhom                % homogenized eigenstress
        timestamp            % timestamp
    end
    properties (Dependent)
        vol                  % volume fractions of RVE phase
        f                    % phase volume fractions related to RVE
    end
    
    methods
        %% CONSTRUCTOR
        function newRVE = RVEData1(namehom,phaseshom,homscheme,shapehom,toli)
            if nargin == 5
                newRVE.name = namehom;
                % validation for phase input
                newRVE.phases=[];
                for i=1:numel(phaseshom)
                    if isa(phaseshom{i},'PhaseData') || isa(phaseshom{i},'PhaseDataE')
                        newRVE.phases.(phaseshom{i}.name) = phaseshom{i};
                    else
                        disp('PHASEHOM{i}')
                        disp(i)
                        disp(phaseshom{i})
                        error('Phase must be a previously defined Phase object, use PhaseData for definition')
                    end
                end
                % validation input for homogenization scheme
                valid_scheme={'MT','SCS','DIF','DIL'};
                if ismember(homscheme,valid_scheme)
                    newRVE.scheme = homscheme;
                end
                
                %newRVE.homorder = homorderhom;
                
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
                %other properties are
            else
                error('6 input properties must be given: name,phases,scheme,shape,tol')
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
        %% CLEAR
        function [obj] = clearhom(obj)
            obj.Chom=[];
            obj.Ainf=[];
            obj.A=[];
            obj.P=[];
            obj.APstroud=[];
            obj.EEinfty=[];
            obj.C0=[];
            obj.timestamp=[];            
        end
        %% ADD/REMOVE/REPLACE PHASES
        function [obj] = addphase(obj,phase_to_add,varargin)
            if ~ismember('noclear',varargin)
            obj=clearhom(obj);
            end
            obj.phases.(phase_to_add.name)=(phase_to_add);        
        end
        function [obj] = rmphase(obj,phase_to_rm,varargin)
            if ~ismember('noclear',varargin)
            obj=clearhom(obj);
            end
            obj.phases=rmfield(obj.phases,phase_to_rm);        
        end
        function [obj] = replacephase(obj,phase_to_rm,phase_to_add,varargin)
            if ~ismember('noclear',varargin)
            obj=clearhom(obj);
            end
            obj.phases=rmfield(obj.phases,phase_to_rm);        
            obj.phases.(phase_to_add.name)=(phase_to_add);        
        end
        %% HOMOGENIZATION
        function [obj,phasehom] = homE(obj,optguess) %% Elastic Homogenization
            if nargin<1.5
                optguess=NaN;
            end
            %gengerate string and function handle for homogenization scheme
            str=['fun_hom_',obj.scheme];
            fhom=str2func(str);
            %evaluate homogenization
            [iC,iA,iAinf,iP,iAPstroud,iEEinfty,iC0]=fhom(obj,optguess);
            obj.Chom=iC;
            obj.Ainf=iAinf;
            obj.A=iA;
            obj.P=iP;
            obj.APstroud=iAPstroud;
            obj.EEinfty=iEEinfty;
            obj.C0=iC0;
            obj.timestamp=now;
            
            % generate new material and phase
            if sum(sum(abs(iC-transpose(iC))))>1e-10
               warning(['homogenized stiffness tensor is not symmetric - something went wrong - sb0=',num2str(sum(sum(abs(iC-transpose(iC)))))])
            end
            Mathom=MatData(obj.name,'hom',(iC+transpose(iC))/2,false,false);
            phasehom=PhaseData(obj.name,Mathom,obj.shape,obj.vol);
        end
        function [obj] = compQ(obj) %% Compute eigenstrain influence tensors
            iQ=fun_eigenstraininfl(obj);
            %NOT TESTED
            obj.Q=iQ;
        end
        function [obj,p] = homPI(obj) %% Homogenization of phase eigenstresses
            sumPIhom=zeros(1,6);
            fld=fieldnames(obj.phases);
            for i=1:numel(fld)
                % check whether all phase eigenstress are defined
                if isempty(obj.phases.(fld{i}).p)
                    warning(['eigenstress of phase ',obj.phases.(fld{i}).name,' not defined - set to zero temporarily'])
                    pip=NaN;
                else
                    % Homogenize eigenstresses
                    pip=obj.phases.(fld{i}).p;                      
                end
                if ~isnan(pip)
                    if ismember(obj.phases.(fld{i}).shape.dis,{'single'})
                        % single inclusion
                        sumPIhom=sumPIhom+obj.f.(fld{i})*pip*obj.A.(fld{i});
                    else
                        % isotropically oriented inclusion
                        % eigenstress need to be isotropically oriented as
                        % well, matrix need to be isotropic %CHECK THIS HERE
                        % check isotropy
                        checkisomat=fun_check_isotropy(obj.C0);
                        if strcmp(checkisomat,'iso')
                            if pip(1)-pip(2)<1e-10 && pip(2)-pip(3)<1e-10
                                sumPIhom=sumPIhom+obj.f.(fld{i})*pip*obj.A.(fld{i});
                            else
                                PIe3=pip*(obj.APstroud.(fld{i}).Ainfe3*obj.EEinfty);
                                PIiso=1/3*(PIe3(1)+PIe3(2)+PIe3(3));
                                sumPIhom=sumPIhom+obj.f.(fld{i})*PIiso;
                            end
                        else
                            error(['Eigenstress of isotropically oriented phase ', obj.phases.(fld{i}).name,' should be isotropic!, Matrix must be isotropic too'])
                        end
                    end
                end
            end
            obj.PIhom=sumPIhom;
            p=obj.PIhom;
        end
        
    end
end