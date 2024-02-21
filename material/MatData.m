classdef MatData
    properties
        name                 % Name of the Material
        iso                  % isotropic, transversally isotropic, orthotropic, anisotropic,
        C                    % Stiffness tensor of the Material [GPa]
        visco                % viscoelastic [1..yes,0..no]
        fail                 % Prone to failure [1..yes,0..no]
    end
    properties (Dependent)
        isoc     
    end
    
    methods
        % CONSTRUCTOR
        function newMat = MatData(material,isiso,stiffnessT,isvisco,doesfail) % constructer
            if nargin > 0
                newMat.name = material;                
                % validation input for isotropy
                valid_iso={'iso','transiso','ortho','aniso','hom','pore'};
                if ismember(isiso,valid_iso)
                    newMat.iso = isiso;
                else
                    error('Isotropy must be one of the following strings: iso,transiso,ortho,aniso');
                end
                
                % validation input for C
                if strcmp(isiso,'iso')&& iscell(stiffnessT) && numel(stiffnessT)==3
                    tmp=strcat('fun_Cfrom',stiffnessT{1});
                    newMat.C=feval(tmp,stiffnessT{2},stiffnessT{3});    
                elseif sum(sum(isnan(stiffnessT)))==0 && sum(size(stiffnessT)==size(magic(6)))==2
                    if issymmetric(stiffnessT)
                        newMat.C = stiffnessT;
                    else
                        error('Stiffness tensor must be symmetric')
                    end
                elseif strcmp(isiso,'hom')
                    newMat.C = zeros(6)+NaN;
                else
                    error('Stiffness Tensor must be 6x6 Matrix, or a 3x1cell if isotropic with syntax {''Enu'',E,nu} or {''kmu'',k,mu}');
                end
                
                % validation for logical inputs visco
                if islogical(isvisco)
                    newMat.visco = isvisco;
                else
                    error('Logical input required for visco')
                end
                
                % validation for logical inputs visco
                valid_failure={'Rankine','Tresca','vonMises','MC','DP'};
                tmp=strcat('Failure is either ''false'', or a cell with structure {''Failcrit'',[par1,par2,...]}, whereby Failcrit is one of the following ',valid_failure);
                if islogical(doesfail);
                    if doesfail==0
                        newMat.fail.i = 0;
                        newMat.fail.crit = '';
                        newMat.fail.par = [];
                    else
                        error(tmp)
                    end
                elseif iscell(doesfail) && numel(doesfail)==2 && sum(sum(strcmp(doesfail{1},valid_failure)))==1
                    newMat.fail.i = 1;
                    newMat.fail.crit = doesfail{1};
                    newMat.fail.par = doesfail{2};
                else
                    failure
                    error(tmp)
                end                
            end
        end
        
        
        % ONLY IF ISOTROPIC
        function isoc = get.isoc(obj)
            if strcmp('iso',obj.iso)==1
                isoc.k=obj.C(1,1)-4/3.*0.5*obj.C(6,6);
                isoc.mu=0.5*obj.C(6,6);
                [isoc.E,isoc.nu]=fun_Enu_from_kmu(isoc.k,isoc.mu);
            end
        end
        
        function obj = set.isoc(obj,~)
            fprintf('%s%d\n','Elastic Properties read as: ',obj.isoc)
            error('You cannot set elastic Properties property directly');
        end
        
        
        % DISPLAY
        function disp(newMat)
            fprintf('******************************\n%stropic Material: %s\n******************************\n',newMat.iso,newMat.name)
            if strcmp(newMat.iso,'iso')
                fprintf('isotropic constants: E=%0.5g, nu=%0.5g, k=%0.5g, mu=%0.5g\n',...
                    newMat.isoc.E,newMat.isoc.nu,newMat.isoc.k,newMat.isoc.mu)
            else
                tmp = [repmat('%6.5g ', 1, size(newMat.C,2)-1), '%6.5g\n'];
                fprintf('Stiffness:\n')
                fprintf(tmp, newMat.C.')
            end
            fprintf('viscoelastic: %i\nfailure: %s\n',newMat.visco,newMat.fail.crit)
            if newMat.fail.i==1
                fprintf('failure criterion: %s, with paramaters: %s',newMat.fail.crit,num2str(newMat.fail.par))
            end
        end
    end
end