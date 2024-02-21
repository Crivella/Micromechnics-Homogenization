function [Ainf,P,addout] = fun_AinfALL_R(phaseobj,C0,checkiso,tol)
% get identity tensors
I4 = fun_I4();

% addition output for isotropically oriented phase
addout.Ainfe3=[];
addout.Pe3=[];
addout.Ainf=[];
addout.P=[];

% R tensor for imperfect IF
if isnan(phaseobj.IF)
    R=zeros(6);
else
    %R=fun_R_ellipsoid(phaseobj.shape.ar,phaseobj.shape.sr,phaseobj.IF(1),phaseobj.IF(2),tol);
    R=fun_R_spheroid(phaseobj.shape.sr,phaseobj.IF(1),phaseobj.IF(2));
end
addout.R=R;

% single inclusion
if ismember(phaseobj.shape.dis,{'single'})
    if ismember(phaseobj.shape.def,{'sphere','matrix'}) || sum(isnan(phaseobj.shape.ori))==2 || (phaseobj.shape.ori(1)==0 && phaseobj.shape.ori(2)==0)
        % evaluate Hill Tensor (in e3-direction)
        P=fun_P(phaseobj,C0,checkiso,tol);
        Ainf=inv(I4.I+P*(phaseobj.mat.C-C0)+(I4.I-P*C0)*R*phaseobj.mat.C);
    else
        % evaluate Hill Tensor (for rotated orientation)
        Q4=fun_Q4_bp(phaseobj.shape.ori(1),phaseobj.shape.ori(2));
        Q4t=transpose(Q4);
        %Q2=fun_Q2_bp(phaseobj.shape.ori(1),phaseobj.shape.ori(2));
        %v=[0,0,1]';
        %v_trans=Q2*v;
        %v_rot=transpose(Q2)*v;
        if strcmp(checkiso,'iso')
            %no stiffness tensor rotation necessary
            P_e3=fun_P(phaseobj,C0,checkiso,tol);
            C0_aniso=C0;
        else
            % Transformation of matrix stiffness components in negative needle orientation
            C0_aniso=Q4*C0*Q4t;
            P_e3=fun_P_ellipsoid_aniso(C0_aniso,phaseobj.shape.ar,phaseobj.shape.sr,tol);
        end
        % Transformation back into needle orientation
        P=(Q4t * P_e3 * Q4);
        R=(Q4t * R * Q4);


        checkisoincl=fun_check_isotropy(phaseobj.mat.C);
        if ~strcmp(checkisoincl,'iso')
            % Transformation of inclusion stiffness components in needle orientation
            C_i=(Q4t * phaseobj.mat.C * Q4);
        else
            C_i=phaseobj.mat.C;
        end
        Ainf=inv(I4.I+P*(C_i-C0)+(I4.I-P*C0)*R*C_i);
    end
end

% isotropically oriented inclusion
if ismember(phaseobj.shape.dis,{'iso'})
    
        %disp(['isotropic averaging for phase ',phaseobj.name])
        % evaluate Hill Tensor (in e3-direction)
        P=fun_P(phaseobj,C0,checkiso,tol);
        Ainf=inv(I4.I+P*(phaseobj.mat.C-C0)+(I4.I-P*C0)*R*phaseobj.mat.C);
        addout.Ainfe3=Ainf;
        addout.Pe3=P;
    
    if ismember(checkiso,{'iso'})
        % A) isotropic averaging
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
            P_e3=fun_P_spheroid_transiso(C0_aniso,phaseobj.shape.ar,phaseobj.shape.sr,tol); %% HERE SHOULD BE ANISO, BUT FOR STRENGTH HOMOGENIZATION ITS SUFFICIENT
            %P_e3=fun_P_ellipsoid_aniso(C0_aniso,phaseobj.shape.ar,phaseobj.shape.sr,tol);
            % Transformation back into needle orientation
            P_i=(Q4t * P_e3 * Q4);
            R_i=(Q4t * R * Q4);
            addout.P(:,:,j)=P_i;

            checkisoincl=fun_check_isotropy(phaseobj.mat.C);
            if ~ismember(checkisoincl,{'iso'})
                % Transformation of inclusion stiffness components in needle orientation
                C_i=(Q4t * phaseobj.mat.C * Q4);
            else
                C_i=phaseobj.mat.C;
            end

            % strain concentration tensors for all stroud directions
            addout.Ainf(:,:,j)=inv(I4.I+P_i*(C_i-C0)+(I4.I-P_i*C0)*R_i*C_i);
            sumAinf=sumAinf+stroud.weight(j)*addout.Ainf(:,:,j);
        end
        Ainf=sumAinf;
    end
end


% n-layered inclusion
if ismember(phaseobj.shape.def,{'nlayered'})
    if ismember(checkiso,{'iso'})
        if ismember(fun_check_isotropy(phaseobj.mat.C),{'iso'})
            if ~exist(addout.nlayered,'var') %setup int_m matrix
                addout.nlayered=zeros(phaseobj.shape{2},3);
            end
            addout.nlayered(phaseobj.shape{2},:)=[phaseobj.mat.isoc.k,phaseobj.mat.isoc.mu,phaseobj.vol];          
        end
    else
        error('n-layered inclusion must be embedded in isotropic matrix')
    end
end

