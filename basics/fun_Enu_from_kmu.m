function [E,nu] = fun_Enu_from_kmu(k,mu)
%calculate Young's modulus (E) and Poisson's ratio (nu), bulk  (k) and
%shear modulus (mu) from stiffness tensor (C)

% transition
nu = (3*k-2*mu)./(6*k+2*mu);
E = 9*k.*mu./(3*k+mu);
end

