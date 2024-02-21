function [Chomi,invEEinftyi] = fun_ChomALL_2(phaseobj,Ainf,APStroudi,checkiso)
Ci=phaseobj.mat.C;
phaseshape=phaseobj.shape;

if ismember(phaseshape.dis,{'single'})
    % single inclusion
    if sum(isnan(phaseshape.ori))==2 || (phaseshape.ori(1)==0 && phaseshape.ori(2)==0) || strcmp(fun_check_isotropy(Ci),'iso')
        % no orientation (inclusion e3-aligned) or isotropic inclusion        
        Chomi=Ci*Ainf;
        invEEinftyi=(eye(6)+APStroudi.R*Ci)*Ainf;
    else
        % rotated (non-isotropic) inclusion
        Q4=fun_Q4_bp(phaseshape.ori(1),phaseshape.ori(2));
        Q4t=transpose(Q4);
        Ci_trans=(Q4t * Ci * Q4);
        R_trans=(Q4t * APStroudi.R * Q4);
        Chomi=Ci_trans*Ainf;
        invEEinftyi=(eye(6)+R_trans*Ci_trans)*Ainf;
    end
else
    % isotropically oriented inclusion
    if ismember(phaseshape.dis,{'iso'})
        if strcmp(fun_check_isotropy(Ci),'iso')
            % isotropic inclusion
            Chomi=Ci*Ainf;
            invEEinftyi=fun_isoav((eye(6)+APStroudi.R*Ci)*APStroudi.Ainfe3);
        elseif  ismember(checkiso,{'iso'})
            % non-isotropic inclusion
            Chomi=fun_isoav(Ci*APStroudi.Ainfe3);
            invEEinftyi=fun_isoav((eye(6)+APStroudi.R*Ci)*APStroudi.Ainfe3);
        else
            % stroud integration, 15 directions
            Chomi=zeros(6);
            stroud=fun_stroud15(); %stroud=fun_stroud28();
            stroud.num=numel(stroud.azi);
            for j=1:stroud.num
                    % Transformation of inclusion stiffness components in needle orientation
                    Q4=fun_Q4_bp(stroud.azi(j),stroud.zeni(j));
                    Q4t=transpose(Q4);
                    Ci_trans=(Q4t * Ci * Q4);
                    R_trans=(Q4t * APStroudi.R * Q4);
                
                % strain concentration tensors for all stroud directions
                Chomi=Chomi+stroud.weight(j)*Ci_trans*APStroudi.Ainf(:,:,j);
                invEEinftyi=invEEinftyi+stroud.weight(j)*(eye(6)+R_trans*Ci_trans)*APStroudi.Ainf(:,:,j);
            end            
        end
    end
end
