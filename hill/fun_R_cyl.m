% function for computation of R-tensor related to "slightly" imperfect
% interfaces for cylindric inclusions in matrix-inclusion problems, 
% see Qu 1993a,b and Dinzart & Sabar 2017
%
%-----------------------------------------------------------------------
% INPUT: 
% alpha ... interface sliding parameter, alpha>=0
% beta ... interface separation parameter, beta>=0
% a ... cylinder radius
%
% OUTPUT: R-tensor
% R ... R-tensor
%
%-----------------------------------------------------------------------
%
function[R]=fun_R_cyl(alpha,beta,a)
if nargin<3
    a=1;
end
alphaA=alpha/a;
betaA=beta/a;

%from Dinzart, I think this is wrong
% R1111=alphaA+3*(betaA-alphaA)/4;
% R1122=(betaA-alphaA)/4;
% R1212=alphaA/2+(betaA-alphaA)/4;
% 
% R=[ R1111,R1122,0,0,0,0;
%     R1122,R1111,0,0,0,0;
%     0,0,0,0,0,0;
%     0,0,0,0,0,0;
%     0,0,0,0,0,0;
%     0,0,0,0,0,2*R1212];

%from Qu, I think this is correct
P1111=3*pi/8;
Q1111=9*pi/32;
P=[         P1111,0,0,0,0,0;...
            0,P1111,0,0,0,0;...
            0,0,0,0,0,0;...
            0,0,0,2/4*P1111,0,0;...
            0,0,0,0,2/4*P1111,0;...
            0,0,0,0,0,P1111];
Q=[         Q1111,Q1111/3,0,0,0,0;...
            Q1111/3,Q1111,0,0,0,0;...
            0,0,0,0,0,0;...
            0,0,0,0,0,0;...
            0,0,0,0,0,0;...
            0,0,0,0,0,2/3*Q1111];
R=alpha*P+(beta-alpha)*Q;
end
