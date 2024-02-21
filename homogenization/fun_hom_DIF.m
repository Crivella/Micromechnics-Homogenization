function [Chom,A,Ainf,P,APstroud,EEinfty,C0] = fun_hom_DIF(obj,~)
%fun_hom_DIF computes homogenized stiffness and eigenstress applying the
%DIFFERENTIAL homogenization scheme

%% list of phases in RVE
fld=fieldnames(obj.phases);

%% get C0 = matrix stiffness
defcell=cell(1,numel(fld));
for i=1:numel(fld)
    defcell{i}=obj.phases.(fld{i}).shape.def;
end
[defmat, posmat] = ismember('matrix',defcell);
if ~defmat
    error('matrix phase in RVE is missing')
end
defcell{posmat}='';
if ismember('matrix',defcell);
    error('only one matrix phase in RVE possible')
end

C0=obj.phases.(fld{posmat}).mat.C;

%% homogenization
disp('Starting homogenization by means of the DIFFERENTIAL SCHEME')
% steps for differential scheme
steps_diff=obj.tol*4;
EEinfty=eye(6);

for j=1:steps_diff
    % check isotropy
    checkiso=fun_check_isotropy(C0); % 1... isotropic, 2... transversaly isotropic in e3, 3... other
    
    % initialization
    Chom=zeros(6);
    % strain concentration tensors
    for i=1:numel(fld)
        % solve Eshelby Problem
        %disp(['Eshelby solution for phase ',obj.phases.(fld{i}).name]);
        [Ainf.(fld{i}),P.(fld{i}),APstroud.(fld{i})]=fun_AinfALL_R(obj.phases.(fld{i}),C0,checkiso,obj.tol);
        %disp(Ainf.(fld{i}))
    end
    
    % homogenized stiffness
    for i=1:numel(fld)
        % RVE-related strain concentration tensor
        A.(fld{i})=Ainf.(fld{i});
        
        % homogenized stiffness
        Chomi = fun_ChomALL(EEinfty,Ainf.(fld{i}),obj.phases.(fld{i}).mat.C,obj.phases.(fld{i}).shape,APstroud.(fld{i}),checkiso);
        Chom=Chom+obj.f.(fld{i})*Chomi;
    end
    C0=Chom;
    disp(['  ---> #',num2str(j),' from ',num2str(steps_diff)])
    %disp(C0)
end
end

