function [E,nu,k,mu] = fun_Enu_from_C(C)
%calculate Young's modulus (E) and Poisson's ratio (nu), bulk  (k) and
%shear modulus (mu) from stiffness tensor (C)

% verify input
if ~isequal(size(C), [6,6])
    error('dimension of matrix must be 6x6')
end

% Lame constants
la=C(1,2);
mu=C(4,4)/2;

% transition
k=la+2/3*mu;
nu = (3*k-2*mu)/(6*k+2*mu);
E = 9*k*mu/(3*k+mu);
end

