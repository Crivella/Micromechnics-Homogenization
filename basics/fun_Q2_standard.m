% function for computation of transformation matrix for transformation of
% 2nd order tensors,
% Definition of base vectors according to Physics
%
%-----------------------------------------------------------------------
% INPUT: azi  ... azimuth angle, reference = e_x
%        zeni ... azimuth angle, reference = e_z
% OUTPUT: Q2-Matrix
%
%-----------------------------------------------------------------------
%
function[Q2_standard]=fun_Q2_standard(azi,zeni)
Q21=cos(azi)*cos(zeni); %Definition Standard!
Q22=sin(azi)*cos(zeni);
Q23=-sin(zeni);
Q31=-sin(azi);
Q32=cos(azi);
Q33=0;
Q11=cos(azi)*sin(zeni);
Q12=sin(azi)*sin(zeni);
Q13=cos(zeni);

Q2_standard=[Q11,Q12,Q13;
             Q21,Q22,Q23;
             Q31,Q32,Q33];
clear Q21 Q22 Q23 Q31 Q32 Q33 Q11 Q12 Q13
end