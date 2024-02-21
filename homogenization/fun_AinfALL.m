function [Ainf,P,addout] = fun_AinfALL(phaseobj,C0,checkiso,tol)
% get identity tensors
I4 = fun_I4();

% addition output for isotropically oriented phase
addout.Ainfe3=[];
addout.Pe3=[];
addout.Ainf=[];
addout.P=[];

% single inclusion
if ismember(phaseobj.shape.dis,{'single'})
    if ismember(phaseobj.shape.def,{'sphere','matrix'}) || sum(isnan(phaseobj.shape.ori))==2 || (phaseobj.shape.ori(1)==0 && phaseobj.shape.ori(2)==0)
        % evaluate Hill Tensor (in e3-direction)
        P=fun_P(phaseobj,C0,checkiso,tol);
        Ainf=inv(I4.I+P*(phaseobj.mat.C-C0));
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

        checkisoincl=fun_check_isotropy(phaseobj.mat.C);
        if ~strcmp(checkisoincl,'iso')
            % Transformation of inclusion stiffness components in needle orientation
            C_i=(Q4t * phaseobj.mat.C * Q4);
        else
            C_i=phaseobj.mat.C;
        end
        Ainf=inv(I4.I+P*(C_i-C0));
    end
end

% isotropically oriented inclusion
if ismember(phaseobj.shape.dis,{'iso'})
    if ismember(checkiso,{'iso'})
        %disp(['isotropic averaging for phase ',phaseobj.name])
        % evaluate Hill Tensor (in e3-direction)
        P=fun_P(phaseobj,C0,checkiso,tol);
        Ainf=inv(I4.I+P*(phaseobj.mat.C-C0));

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
                P=NaN;
                Ainf=NaN: % calculate after the setup of the in_m matrix
                %in_m=[  k1, mu1, R1;... %iso Fe core
                %    k2, mu2, R2;... %iso FeCO3 shell
                %    k0, mu0, 100];
                %sol_int=fun_HZ_int(in_m);
                %Ainf=sol_int(1,5)*I4.vol+sol_int(1,7)*I4.dev;
                %Ainf=sol_int(2,5)*I4.vol+sol_int(2,7)*I4.dev;
            end
        end
    else
        error('n-layered inclusion must be embedded in isotropic matrix')
    end

end

