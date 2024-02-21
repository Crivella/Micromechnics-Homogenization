% function for computation of R-tensor related to "slightly" imperfect
% interfaces for cylindric inclusions in matrix-inclusion problems, 
% see Qu 1993a,b and Dinzart & Sabar 2017
%
%-----------------------------------------------------------------------
% INPUT: 
% alpha ... interface sliding parameter, alpha>=0
% beta ... interface separation parameter, beta>=0
% a ... sphere radius
%
% OUTPUT: R-tensor
% R ... R-tensor
%
%-----------------------------------------------------------------------
%
function[R]=fun_R_sph(alpha,beta,a)
if nargin<3
    a=1;
end
alphaA=alpha/a;
betaA=beta/a;
Rvol=betaA;
Rdev=(3*alphaA+2*betaA)/5;

I=[1 0 0 0 0 0;...
   0 1 0 0 0 0;...
   0 0 1 0 0 0;...
   0 0 0 1 0 0;...
   0 0 0 0 1 0;...
   0 0 0 0 0 1];

K=[1/3 1/3 1/3 0 0 0;...
   1/3 1/3 1/3 0 0 0;...
   1/3 1/3 1/3 0 0 0;...
   0   0   0   0 0 0;...
   0   0   0   0 0 0;...
   0   0   0   0 0 0];

J=I-K;

R=Rvol*K+Rdev*J;
end
