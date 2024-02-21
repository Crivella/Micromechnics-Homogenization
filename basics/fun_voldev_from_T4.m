function [Avol,Adev] = fun_voldev_from_T4(T4)
%calculate volumetric and deviatoric part of foruth order tensor
% verify input
if ~isequal(size(T4), [6,6])
    error('dimension of matrix must be 6x6')
end

% Lame constants
la=T4(1,2);
mu=T4(4,4)/2;

% transition
Adev=2*mu;
Avol=3*la+Adev;
end

