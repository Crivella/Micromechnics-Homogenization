classdef PhaseData < matlab.mixin.Copyable
    properties
        name                 % Name of the Phase
        mat                  % Material Object
        shape                % Shape
        vol                  % Phase volume fraction related to macromaterial
        p                    % eigenstress
        IF                   % interface compliance
    end
    properties (Dependent)
        %orientation          % Phase orientation
    end
    events
        %defhomMAT
    end
    
    methods
        % CONSTRUCTOR
        function newphase = PhaseData(phasename,matobj,phaseshape,phasevol,phasepi,IF)
            if nargin > 0
                newphase.name = phasename;
                % validation input for matobj
                if isa(matobj,'MatData')
                    newphase.mat = matobj;
                else
                    error('Material must be a previously defined Material object, use MatData for definition')
                end
                % validation input for shape
                valid_shape={'matrix','sphere','spheroid','ellipsoid','nlayered'};
                if sum(sum(strcmp(phaseshape{1},valid_shape)))==1
                    newphase.shape.def = phaseshape{1};
                    newphase.shape.ar = 1;
                    newphase.shape.sr = 1;
                    newphase.shape.dis = 'single';
                    newphase.shape.ori = [NaN,NaN];
                    newphase.shape.layer = 0;
                else
                    error(strjoin({'shape must be a cell array, whereby the first element is either of the following strings',strjoin(valid_shape)}));
                end
                if sum(sum(strcmp(phaseshape{1},{'sphere'})))==1 && numel(phaseshape)==2
                    if sum(sum(strcmp(phaseshape{2},'iso')))==1
                        newphase.shape.dis = 'iso';
                    elseif  numel(phaseshape{2})==2 && sum(isnan(phaseshape{3}))==0;
                        newphase.shape.ori = phaseshape{3};
                    else
                        error('orientation must be either the string iso; or a 2-element vector with azimuth, zenith angles');
                    end
                end
                if sum(sum(strcmp(phaseshape{1},{'sphere','matrix','nlayered'})))~=1
                    % ONLY IF NONSPHERICAL
                    if sum(sum(strcmp(phaseshape{1},{'spheroid'})))==1
                        newphase.shape.sr=phaseshape{2};
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
                if sum(sum(strcmp(phaseshape{1},{'nlayered'})))==1
                    newphase.shape.layer=phaseshape{2};
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
            end
            % eigenstress
            if nargin<5
                newphase.p=NaN;
                newphase.IF=NaN;
            else
                %Input validation for eigenstress pi
                if size(phasepi)==[1,1]
                    if phasepi==0
                        newphase.p=NaN;
                    else
                        newphase.p=fun_std2comp(eye(3)*phasepi); %isotropic stress
                    end
                elseif size(phasepi)==[3,3]
                    newphase.p=fun_std2comp(phasepi);
                elseif size(phasepi)==[1,6]
                    newphase.p=phasepi;
                elseif size(phasepi)==[6,1]
                    newphase.p=phasepi';
                else
                    error('input must be either scalar (isotropic eigenstress), or tensor in 3x§ or 6x1 notation')
                end
            end
            % eigenstress
            if nargin<6
                newphase.IF=NaN;
            else
                newphase.IF=IF;
            end
        end
        
    end 
end