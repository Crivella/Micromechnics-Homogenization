% C = fun_CfromEnu(E,nu) returns the isotropic elastic stiffness tensor of 
% a material with Young's modulus E and Poisson's ratio nu
%-----------------------------------------------------------------------
% INPUT: E  ... Young's modulus ->scalar
%        nu ... Poisson's ratio ->scalar
% OUTPUT: C ... stiffness tensor ->4th tensor
%-----------------------------------------------------------------------
%
function[C]=fun_CfromEnu(E,nu)
[k,mu]=fun_kmu_from_Enu(E,nu);
C=fun_Cfromkmu(k,mu);
end
