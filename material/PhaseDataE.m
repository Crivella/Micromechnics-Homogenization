classdef PhaseDataE < handle
    properties
        name                 % Name of the Phase
        mat                  % Material Object
        shape                % Shape
        vol                  % Phase volume fraction related to macromaterial
        pi                   % Phase eigenstress (global base)
    end
    properties (Dependent)
        %orientation          % Phase orientation
    end
    events
        %defhomMAT
    end
    
    methods
        % CONSTRUCTOR
        function newphase = PhaseDataE(phasename,matobj,phaseshape,phasevol,phasepi)
            if nargin > 0
                newphase.name = phasename;
                % validation input for matobj
                if isa(matobj,'MatData')
                    newphase.mat = matobj;
                else
                    error('Material must be a previously defined Material object, use MatData for definition')
                end
                % validation input for shape
                valid_shape={'matrix','sphere','needle','disc','spheroid','ellipsoid'};
                if sum(sum(strcmp(phaseshape{1},valid_shape)))==1
                    newphase.shape.def = phaseshape{1};
                    newphase.shape.ar = 1;
                    newphase.shape.sr = 1;
                    newphase.shape.dis = 'single';
                    newphase.shape.ori = [NaN,NaN];
                else
                    error(strjoin({'shape must be a cell array, whereby the first element is either of the following strings',strjoin(valid_shape)}));
                end
                if sum(sum(strcmp(phaseshape{1},{'sphere','matrix'})))~=1
                    % ONLY IF NONSPHERICAL
                    if sum(sum(strcmp(phaseshape{1},{'spheroid'})))==1
                        newphase.shape.ar=phaseshape{2};
                    end
                    if sum(sum(strcmp(phaseshape{1},{'ellipsoid'})))==1
                        newphase.shape.ar=phaseshape{2}(1);
                        newphase.shape.sr=phaseshape{2}(2);
                    end
                    % orientation
                    if numel(phaseshape)==3
                    if sum(sum(strcmp(phaseshape{3},'iso')))==1
                        newphase.shape.dis = 'iso';
                    elseif  numel(phaseshape{3})==2 && sum(isnan(phaseshape{3}))==0;
                        newphase.shape.ori = phaseshape{3};
                    else
                        error('orientation must be either the string iso; or a 2-element vector with azimuth, zenith angles');
                    end
                    end
                end
                %newphase.RVE = RVEname;
                %                   if ~isobject(RVEname)
                %                       notify(obj,'defhomMAT')
                %                       disp('LALALA')
                %                   end
                
                if numel(phasevol)==1 && phasevol<=1+1e-6
                    newphase.vol = phasevol;
                else
                    error('volume fraction must be a number smaller or equal than 1')
                end
                
                %Input validation for eigenstress pi
                if size(phasepi)==[1,1]
                    if phasepi==0
                        newphase.pi=NaN;
                    else
                        newphase.pi=fun_std2comp(eye(3)*phasepi); %isotropic stress
                    end
                elseif size(phasepi)==[3,3]
                    newphase.pi=fun_std2comp(phasepi);
                elseif size(phasepi)==[1,6]
                    newphase.pi=phasepi;
                elseif size(phasepi)==[6,1]
                    newphase.pi=phasepi';
                else
                    error('input must be either scalar (isotropic eigenstress), or tensor in 3x§ or 6x1 notation')
                end
                
                    
            end
        end
        
        %         % ONLY IF NONSPHERICAL
        %         function orientation = get.isoc(obj)
        %             if strcmp('iso',obj.iso)==1
        %                 isoc.k=obj.C(1,1)-4/3.*0.5*obj.C(6,6);
        %                 isoc.mu=0.5*obj.C(6,6);
        %                 [isoc.E,isoc.nu]=fun_Enu_from_kmu(isoc.k,isoc.mu);
        %             end
        %         end
        %
        %         function obj = set.isoc(obj,~)
        %             fprintf('%s%d\n','Elastic Properties read as: ',obj.isoc)
        %             error('You cannot set elastic Properties property directly');
        %         end
    end
    
    
    
    
end