% function for computation of transformation matrix for transformation of
% 2th order tensors,
% Definition of base vectors according to Bernhard Pichler, thus for
% azi=zeni=0, the e_x,e_y,e_z base vectors are equal to the base vectors
% e_zeni, e_azi, e_r
% This defnition deviates from standard definition, but allows to directly
% use this matrix for stroud integration
%
%-----------------------------------------------------------------------
% INPUT: azi  ... azimuth angle, reference = e_x
%        zeni ... azimuth angle, reference = e_z
% OUTPUT: Q2-Matrix
%
%-----------------------------------------------------------------------
%
function[Q2_bp]=fun_Q2_bp(azi,zeni)
Q11=cos(azi)*cos(zeni); %Definition Bernhard!
Q12=sin(azi)*cos(zeni);
Q13=-sin(zeni);
Q21=-sin(azi);
Q22=cos(azi);
Q23=0;
Q31=cos(azi)*sin(zeni);
Q32=sin(azi)*sin(zeni);
Q33=cos(zeni);

Q2_bp=[ Q11,Q12,Q13;...
        Q21,Q22,Q23;...
        Q31,Q32,Q33];
clear Q21 Q22 Q23 Q31 Q32 Q33 Q11 Q12 Q13
end