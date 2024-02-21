function [Ainf,P,addout] = fun_AinfALL2(phaseobj,C0,checkiso,tol)
% get identity tensors
I4 = fun_I4();

% addition output for isotropically oriented phase
addout.Ainfe3=[];
addout.Pe3=[];
addout.Ainf=[];
addout.P=[];

% single inclusion or isotropically oriented inclusion in isotropic matrix
if ismember(phaseobj.shape.dis,{'single'}) || (ismember(phaseobj.shape.dis,{'iso'}) && ismember(checkiso,{'iso'}))
    if ismember(phaseobj.shape.def,{'sphere','matrix'}) || sum(isnan(phaseobj.shape.ori))==2 || (phaseobj.shape.ori(1)==0 && phaseobj.shape.ori(2)==0) || ismember(phaseobj.shape.dis,{'iso'})
        % evaluate Hill Tensor (in e3-direction)
        P=fun_P(phaseobj,C0,checkiso,tol);
    else
        % evaluate Hill Tensor (for rotated orientation)
        Q4=fun_Q4_bp(phaseobj.shape.ori(1),phaseobj.shape.ori(2));
        Q4t=transpose(Q4);
        if strcmp(checkiso,'iso')
            %no stiffness tensor rotation necessary
            P_e3=fun_P_OLD(phaseobj,C0,checkiso,tol);
        else
            % Transformation of matrix stiffness components in negative needle orientation
            C0_aniso=Q4*C0*Q4t;
            P_e3=fun_P_ellipsoid_aniso(C0_aniso,phaseobj.shape.ar,phaseobj.shape.sr,tol);
        end
        % Transformation back into needle orientation
        P=(Q4t * P_e3 * Q4);
    end
    
    % Eshelby-Problem-related strain concentration tensor
    Ainf=inv(I4.I+P_e3*(phaseobj.mat.C-C0));
end

% isotropically oriented inclusion
if ismember(phaseobj.shape.dis,{'iso'})
    if ismember(checkiso,{'iso'})
        %disp(['isotropic averaging for phase ',phaseobj.name])
        % A) isotropic averaging
        addout.Ainfe3=Ainf;
        addout.Pe3=P;
        Ainf=fun_isoav(Ainf); % isotropic averaging
        
    else

        %disp(['Stroud integration for phase ',phaseobj.name])
        % B) stroud integration
        if tol<9
            stroud=fun_stroud15(); 
        elseif tol<12
            stroud=fun_stroud28();
        elseif tol==12
            stroud=getLebedevHemiSphere(86);
        elseif tol==13
            stroud=getLebedevHemiSphere(146);
        elseif tol==14
            stroud=getLebedevHemiSphere(350);
        elseif tol==15
            stroud=getLebedevHemiSphere(1202);
        elseif tol==16
            stroud=getLebedevHemiSphere(5810);
        end
        stroud.num=numel(stroud.azi);
        sumAinf=zeros(6);
        addout.Ainf=zeros(6,6,stroud.num);
        addout.P=zeros(6,6,stroud.num);
        P=zeros(6);
        
        for j=1:stroud.num
            % Transformation matrices
            Q4=fun_Q4_bp(stroud.azi(j),stroud.zeni(j));
            Q4t=transpose(Q4);
            
            % Transformation of matrix stiffness components in negative needle orientation
            C0_aniso=Q4*C0*Q4t;
            % Hill tensor for that orienation
            P_e3=fun_P_spheroid_transiso(C0_aniso,phaseobj.shape.ar,phaseobj.shape.sr,tol); %% HERE SHOULD BE ANISO, BUT FOR STRENGTH HOMOGENIZATION IT WORKS
            %P_e3=fun_P_ellipsoid_aniso(C0_aniso,phaseobj.shape.ar,phaseobj.shape.sr,tol);
            % Transformation back into needle orientation
            P_i=(Q4t * P_e3 * Q4);
            addout.P(:,:,j)=P_i;
            
            checkisoincl=fun_check_isotropy(phaseobj.mat.C);
            if ~ismember(checkisoincl,{'iso'})
                % Transformation of inclusion stiffness components in needle orientation
                C_i=(Q4t * phaseobj.mat.C * Q4);
            else
                C_i=phaseobj.mat.C;
            end
            
            % strain concentration tensors for all stroud directions
            addout.Ainf(:,:,j)=inv(I4.I+P_i*(C_i-C0));
            sumAinf=sumAinf+stroud.weight(j)*addout.Ainf(:,:,j);
        end
        Ainf=sumAinf;
    end
end

% n-layered inclusion
if ismember(phaseobj.shape.def,{'nlayered'})
    if ismember(checkiso,{'iso'})
        if ismember(fun_check_isotropy(phaseobj.mat.C),{'iso'})
            if phaseobj.shape.layer.n==1
                addout.nlayered=1; %calculate the integration coefficients for first layer
            end
        end
    else
        error('n-layered inclusion must be embedded in isotropic matrix')
    end

end

