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

%%%%% START OF GENERAL HOMOGENIZATION PART-> TO BE FUNCTIONALIZED %%%%%%%%%
% initialization
invEEinfty=zeros(6);
Chom=zeros(6);
for i=1:numel(fld)
    % Eshelby Problem: (infinite) strain concentration tensors
    [Ainf.(fld{i}),P.(fld{i}),APstroud.(fld{i})]=fun_AinfALL_R(obj.phases.(fld{i}),C0,checkiso,obj.tol);
    % % RVE: homogenized stiffness
    [Chomi,invEEinftyi] = fun_ChomALL_2(obj.phases.(fld{i}),Ainf.(fld{i}),APstroud.(fld{i}),checkiso);
    Chom=Chom+obj.f.(fld{i})*Chomi;
    invEEinfty=invEEinfty+obj.f.(fld{i})*invEEinftyi;
end
Chom=Chom/invEEinfty;
for i=1:numel(fld)
    % RVE-related strain concentration tensor
    A.(fld{i})=Ainf.(fld{i})/invEEinfty;
end
EEinfty=inv(invEEinfty);
%%%%% END OF GENERAL HOMOGENIZATION PART-> TO BE FUNCTIONALIZED %%%%%%%%%



