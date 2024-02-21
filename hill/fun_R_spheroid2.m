% function for computation of R-tensor related to "slightly" imperfect
% interfaces for spheroidal and prolate inclusions in matrix-inclusion problems,
% see Lee 2018
% take care: in the paper x1 is the axis of symmetry, here x3, permutation
% necessary
%-----------------------------------------------------------------------
% INPUT:
% asp   . aspect ratio (a1./a2)
% slend . slenderness ratio (a1./a3)
% alpha . interface sliding parameter, alpha>=0 [m/GPa]
% beta . interface separation parameter, beta>=0 [m/GPa]
% tol . tolerance integer

%
% OUTPUT: R-tensor
% R .... R-tensor
%
%-----------------------------------------------------------------------
%
function[R]=fun_R_spheroid2(sr,alpha,beta)
rho=1/sr; %aspect ratio
c=1;%5*1e-6;

if rho>1 %prolate
    P1111=3/4/c*(rho/sqrt(rho^2-1)*(rho^2-2)/(2*(rho^2-1))*asin(sqrt(rho^2-1)/rho)+rho/(2*(rho^2-1)));
    P3333=3/2/c*(rho/(rho^2-1)^(3/2))*asin(sqrt(rho^2-1)/rho)-1./(rho*(rho^2-1));
    P1313=1/4*(P1111+P3333);
    P1212=1/2*P1111;

    Q3333=3/2/c*((2*rho^2+1)/(rho*(rho^2-1)^2)-3*rho/(rho^2-1)^(5/2)*asin(sqrt(rho^2-1)/rho));
    Q1133=3/4/c*(rho*(rho^2+2)/(rho^2-1)^(5/2)*asin(sqrt(rho^2-1)/rho)-3*rho/(rho^2-1)^2);
    Q1111=3/2/c*(rho*(2+rho^2)/(rho^2-1)^2+rho^3*(rho^2-4)/(rho^2-1)^(5/2)*asin(sqrt(rho^2-1)/rho));    
    Q1122=1/3*Q1111;
    P=[P1111,0,0,0,0,0;
        0,P1111,0,0,0,0;
        0,0,P3333,0,0,0;
        0,0,0,2.*P1313,0,0;
        0,0,0,0,2.*P1313,0;
        0,0,0,0,0,2.*P1212];
    Q=[Q1111,Q1122,Q1133,0,0,0;
        Q1122,Q1111,Q1133,0,0,0;
        Q1133,Q1133,Q3333,0,0,0;
        0,0,0,0,0,0;
        0,0,0,0,0,0;
        0,0,0,0,0,0];
    R=alpha*P+(beta-alpha)*Q;
elseif rho<1 %oblate
    
elseif rho==1
    alphaA=alpha/c;
    betaA=beta/c;
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
else
    error()
end
