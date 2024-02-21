% function for computation of R-tensor related to "slightly" imperfect
% interfaces for spheroidal and prolate inclusions in matrix-inclusion problems,
% see Rao 2017 Appendix B
%
%-----------------------------------------------------------------------
% NOTE: spheroid with half axis length a1=a2, a3 (symmetry along x1-x2)
% INPUT:
% slend .... slenderness ratio (a1./a3), >1 oblate, <1 prolate
% alpha .... interface sliding parameter, alpha_a>=0 [1/GPa]
% beta .... interface separation parameter, beta_a>=0 [1/GPa]
% ATTENTION: .... interface sliding parameter need to be normalized with
% the half (!) axis length in the symmetry plane , so alpha and beta are
% actually alpha/a1 and beta/a1

function [R] = fun_R_spheroid(sr,alpha,beta)
rho=1/sr;
if rho==1
    Rvol=alpha;
    Rdev=(3*alpha+2*beta)/5;

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
    if rho>1 %prolate
        P3333=3/2*(rho/((rho^2-1)^(3/2))*asin(sqrt(rho^2-1)/rho)-1/(rho*(rho^2-1)));
        P1111=3/4*(rho*(rho^2-2)/((rho^2-1)^(3/2))*asin(sqrt(rho^2-1)/rho)+rho/(rho^2-1));
        Q3333=3/2*((2*rho^2+1)/(rho*(rho^2-1)^2)-3*rho/((rho^2-1)^(5/2))*asin(sqrt(rho^2-1)/rho));
        Q1133=3/4*(rho*(rho^2+2)/((rho^2-1)^(5/2))*asin(sqrt(rho^2-1)/rho)-3*rho/((rho^2-1)^2));
        Q1111=9/16*(rho*(2+rho^2)/((rho^2-1)^2)+rho^3*(rho^2-4)/((rho^2-1)^(5/2))*asin(sqrt(rho^2-1)/rho));
    elseif rho<1 %oblate
        P3333=3/2*(rho/((rho^2-1)^(3/2))*asinh(sqrt(rho^2-1)/rho)-1/(rho*(rho^2-1)));
        P1111=3/4*(rho*(rho^2-2)/((rho^2-1)^(3/2))*asinh(sqrt(rho^2-1)/rho)+rho/(rho^2-1));
        Q3333=3/2*((2*rho^2+1)/(rho*(rho^2-1)^2)-3*rho/((rho^2-1)^(5/2))*asinh(sqrt(rho^2-1)/rho));
        Q1133=3/4*(rho*(rho^2+2)/((rho^2-1)^(5/2))*asinh(sqrt(rho^2-1)/rho)-3*rho/((rho^2-1)^2));
        Q1111=9/16*(rho*(2+rho^2)/((rho^2-1)^2)+rho^3*(rho^2-4)/((rho^2-1)^(5/2))*asinh(sqrt(rho^2-1)/rho));
    else
        error()
    end
    P1313=1/4*(P3333);
    P1212=1/2*P1111;
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
end
end
