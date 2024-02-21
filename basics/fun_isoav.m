% function for computation of the isotropic average of a 4th order tensor
%
%-----------------------------------------------------------------------
% INPUT: tensT .... fourth-order tensor with minor symmetries in Kelvin-Mandel notation
% OUTPUT: isotropic average of tensT in Kelvin Mandel notation
%
%-----------------------------------------------------------------------
%
function[isoav]=fun_isoav(tensT)
% get unity tensors
I4=fun_I4();

% Isotropic averaging, see Sadowski et al 2015
Fiijj=  tensT(1,1)+tensT(1,2)+tensT(1,3)+...
        tensT(2,1)+tensT(2,2)+tensT(2,3)+...
        tensT(3,1)+tensT(3,2)+tensT(3,3);
    
Fijij=  tensT(1,1)+1/2*tensT(6,6)+1/2*tensT(5,5)+...
        1/2*tensT(6,6)+tensT(2,2)+1/2*tensT(4,4)+...
        1/2*tensT(5,5)+1/2*tensT(4,4)+tensT(3,3);
    
isoav=1/3*Fiijj*I4.vol+1/5*(Fijij-1/3*Fiijj)*I4.dev;
end
