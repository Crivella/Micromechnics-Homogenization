% function for computation of eigenvalues of second order tensors in
% Kelvin-Mandel notation
%
%-----------------------------------------------------------------------
% INPUT: tens2_KV  ... symm. second order tensor in Kelvin Mandel Notation
% OUTPUT: eig_tens2_KV ... eigenvalues of the tensor
%-----------------------------------------------------------------------
%
function[eig_tens2_KV]=fun_eigKV(tens2_KV)

tens2_matrix=[tens2_KV(1),tens2_KV(6),tens2_KV(5);tens2_KV(6),tens2_KV(2),tens2_KV(4);tens2_KV(5),tens2_KV(4),tens2_KV(3)];
eig_tens2_KV=eig(tens2_matrix);