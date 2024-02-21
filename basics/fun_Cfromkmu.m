% C = fun_CfromEnu(E,nu) returns the isotropic elastic stiffness tensor of 
% a material with bulk modulus k and shear modulus mu
%-----------------------------------------------------------------------
% INPUT: k  ... bulk modulus ->scalar
%        mu ... Shear modulus ->scalar
% OUTPUT: C ... stiffness tensor ->4th tensor
%-----------------------------------------------------------------------
%
function[C]=fun_CfromEnu(k,mu)
I4=fun_I4();
C=3*k*I4.vol+2*mu*I4.dev;
end
