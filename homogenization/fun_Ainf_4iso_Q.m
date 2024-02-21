function [RVEobj] = fun_Ainf_4iso_Q(phaseobj,RVEobj)
% NOT TESTED AT ALL

% get identity tensors
I4 = fun_I4();

% check isotropy
checkisomat=fun_check_isotropy(RVEobj.C0);
checkisoinc=fun_check_isotropy(phaseobj.mat.C);

% isotropically oriented inclusion
if ismember(phaseobj.shape.dis,{'iso'})
    % stroud integration
    stroud=fun_stroud15(); %stroud=fun_stroud28();
    stroud.num=numel(stroud.azi);
    addout.Ainf=zeros(6,6,stroud.num);
    addout.P=zeros(6,6,stroud.num);
    for j=1:stroud.num
        % Transformation matrices
        Q4=fun_Q4_bp(stroud.azi(j),stroud.zeni(j));
        Q4t=transpose(Q4);
        
        if strcmp(checkisomat,'iso')
            %no stiffness tensor rotation necessary
            P_e3=fun_P(phaseobj,RVEobj.C0,checkisomat,tol);
        else
            % Transformation of matrix stiffness components in negative needle orientation
            C0_aniso=Q4*C0*Q4t;
            % Hill tensor for that orienation
            P_e3=fun_P_ellipsoid_aniso(C0_aniso,phaseobj.shape.ar,phaseobj.shape.sr,tol);
            % Transformation back into needle orientation
            P_i=(Q4t * P_e3 * Q4);
            addout.P(:,:,j)=P_i;
            
            if ~ismember(checkisoinc,{'iso','transiso'})
                % Transformation of inclusion stiffness components in needle orientation
                C_i=(Q4t * phaseobj.mat.C * Q4);
            else
                C_i=phaseobj.mat.C;
            end
            
            % strain concentration tensors for all stroud directions
            addout.Ainf(:,:,j)=inv(I4.I+P_i*(C_i-C0));
        end
    end
end
RVEobj.APstroud.(phaseobj.name).Ainf=addout.Ainf;
RVEobj.APstroud.(phaseobj.name).P=addout.P;
end

