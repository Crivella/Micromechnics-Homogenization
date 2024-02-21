function [Chomi] = fun_ChomALL(EEinfty,Ainf,Ci,phaseshape,APStroudi,checkiso)

if ismember(phaseshape.dis,{'single'})
    % single inclusion
    if sum(isnan(phaseshape.ori))==2 || (phaseshape.ori(1)==0 && phaseshape.ori(2)==0) || strcmp(fun_check_isotropy(Ci),'iso')
        % no orientation (inclusion e3-aligned) or isotropic inclusion
        Chomi=(Ci*Ainf)*EEinfty;
    else
        % rotated inclusion
        Q4=fun_Q4_bp(phaseshape.ori(1),phaseshape.ori(2));
        Q4t=transpose(Q4);
        Ci_trans=(Q4t * Ci * Q4);
        Chomi=(Ci_trans)*(Ainf*EEinfty);
    end
else
    % isotropically oriented inclusion
    if ismember(phaseshape.dis,{'iso'})
        if strcmp(fun_check_isotropy(Ci),'iso')
            % isotropic inclusion
            %disp('gemma')
            Chomi=(Ci*Ainf)*EEinfty;
        elseif  ismember(checkiso,{'iso'})
            % non-isotropic inclusion
            Chomi=fun_isoav(Ci*APStroudi.Ainfe3)*EEinfty;
        else
            % stroud integration, 15 directions
            Chomi=zeros(6);
            stroud=fun_stroud15(); %stroud=fun_stroud28();
            stroud.num=numel(stroud.azi);
            for j=1:stroud.num
                if ~ismember(fun_check_isotropy(Ci),{'iso'})
                    % Transformation of inclusion stiffness components in needle orientation
                    Q4=fun_Q4_bp(stroud.azi(j),stroud.zeni(j));
                    Q4t=transpose(Q4);
                    Ci_trans=(Q4t * Ci * Q4);
                else
                    Ci_trans=Ci;
                end
                
                % strain concentration tensors for all stroud directions
                Chomi=Chomi+((stroud.weight(j)*Ci_trans*APStroudi.Ainf(:,:,j))*EEinfty);
            end
            
        end
    end
end
