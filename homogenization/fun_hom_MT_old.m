function [Chom,A,Ainf,P,APstroud,EEinfty,C0] = fun_hom_MT(obj,~)
%fun_hom_MT computes homogenized stiffness and eigenstress applying the
%Mori-Tanaka homogenization scheme

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
%disp(['Starting homogenization of ',obj.name,' by means of the MORI-TANAKA SCHEME'])
% check isotropy
checkiso=fun_check_isotropy(C0);

% initialization
invEEinfty=zeros(6);
Chom=zeros(6);

% strain concentration tensors
for i=1:numel(fld)
    % solve Eshelby Problem
    %disp(['Eshelby solution for phase ',obj.phases.(fld{i}).name]);
    [Ainf.(fld{i}),P.(fld{i}),APstroud.(fld{i})]=fun_AinfALL_R(obj.phases.(fld{i}),C0,checkiso,obj.tol);
    
    % link to RVE
    % R tensor for imperfect IF
    if isnan(obj.phases.(fld{i}).IF)
        invEEinfty=invEEinfty+obj.f.(fld{i})*Ainf.(fld{i});
    else
        if ismember(phaseshape.dis,{'iso'})
            invEEinfty=invEEinfty+obj.f.(fld{i})*Ainf.(fld{i})+obj.f.(fld{i})*fun_isoav(APstroud.(fld{i}).R*obj.phases.(fld{i}).mat.C*APstroud.(fld{i}).Ainfe3);
        end
    end
end
%n-layered

EEinfty=inv(invEEinfty);

% homogenized stiffness
for i=1:numel(fld)
    % RVE-related strain concentration tensor
    A.(fld{i})=Ainf.(fld{i})/invEEinfty;
    
    % homogenized stiffness
    Chomi = fun_ChomALL(EEinfty,Ainf.(fld{i}),obj.phases.(fld{i}).mat.C,obj.phases.(fld{i}).shape,APstroud.(fld{i}),checkiso);
    Chom=Chom+obj.f.(fld{i})*Chomi;
end

end

