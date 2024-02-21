% function for computation of the Kronecker Delta
%
%-----------------------------------------------------------------------
% INPUT: i,j .... number
% OUTPUT: 0 or 1 (the latter only if i=j)
%
%-----------------------------------------------------------------------
%
function[kroneckerDelta]=fun_kronecker(i,j)
if i==j
    kroneckerDelta=1;
else
    kroneckerDelta=0;
end
end