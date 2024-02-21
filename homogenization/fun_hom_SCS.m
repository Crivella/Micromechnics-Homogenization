function [Chom,A,Ainf,P,APstroud,EEinfty,C0] = fun_hom_SCS(obj,optguess)
%fun_hom_SCS computes homogenized stiffness and eigenstress applying the
%self-consistent homogenization scheme

%% list of phases in RVE
fld=fieldnames(obj.phases);

%% check if there is a matrix phase
defcell=cell(1,numel(fld));
for i=1:numel(fld)
    defcell{i}=obj.phases.(fld{i}).shape.def;
end
[defmat, ~] = ismember('matrix',defcell);
if defmat
    error('there is a matrix phase in the RVE')
end

%% Precision
mintol=-1;
maxtol=-16;
tol_list=logspace(mintol,maxtol,16);
tolerance=tol_list(obj.tol)/100;

%% initial guess on homogenized stiffness based on rule of mixture

C0=zeros(6);
checkC0=cell(numel(fld),1);

for i=1:numel(fld)
    checkiso=fun_check_isotropy(obj.phases.(fld{i}).mat.C);
    if (strcmp(checkiso,'iso') && strcmp(obj.phases.(fld{i}).shape.def,'sphere')) || strcmp(obj.phases.(fld{i}).shape.dis,'iso')
        % isotropic contribution of phase i
        C0i=obj.f.(fld{i})*fun_isoav(obj.phases.(fld{i}).mat.C);
        checkC0{i}='iso';
    elseif (ismember(checkiso,{'iso','transiso'}) && ismember(obj.phases.(fld{i}).shape.def,{'spheroid','sphere'})) && (sum(isnan(obj.phases.(fld{i}).shape.ori))==2 || (obj.phases.(fld{i}).shape.ori(1)==0 && obj.phases.(fld{i}).shape.ori(2)==0))
        % transversal isotropic contribution of phase i
        C0i=obj.f.(fld{i})*obj.phases.(fld{i}).mat.C;
        checkC0{i}='transiso';
    else
        % anisotropic contribution of phase i
        C0i=obj.f.(fld{i})*obj.phases.(fld{i}).mat.C;
        checkC0{i}='aniso';
    end
    C0=C0+C0i;
end
if sum(isnan(optguess))==0
    C0=optguess;
end

% type of anisotropy of C0
checkiso='iso';
if ismember('transiso',checkC0);
    checkiso='transiso';
elseif ismember('aniso',checkC0);
    checkiso='aniso';
end
%disp(checkiso)
%% homogenization
sumerr=1;
maxcounter=obj.tol^1.5;
%disp(['Starting homogenization of ',obj.name,' by means of the SELF-CONSISTENT SCHEME'])
counter=0;

while sumerr > tolerance && counter<maxcounter
    counter=counter+1;

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

    sumerr=sum(sum(abs(Chom-C0)))/mean(mean(Chom(1:3,1:3)));
    C0=Chom;
    %make C0 transiso
    %C0=(C0+transpose(C0))/2;
    %disp(['  ---> #',num2str(counter),', error = ',num2str(sumerr),'>',num2str(tolerance)])
end
C0=Chom;
end

