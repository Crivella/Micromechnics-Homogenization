% function for computation of P-tensor for spheroidal isotropic inclusions
% aligned with the e3-axis, embedded in isotropic matrix material
%
%-----------------------------------------------------------------------
% INPUT: matrix stiffness tensor C0 (isotropic)
%        slenderness ratio
%
% OUTPUT: P-tensor
%
%-----------------------------------------------------------------------
%
function[P]=fun_P_spheroid_iso(C0,~,sr)
sr=1/sr;
if sr<1 && sr>0 % prolate
    g=sr/(sr^2-1)^(3/2)*(sr*(sr^2-1)^(1/2)-acosh(sr));
elseif sr>1 % oblate   
    g=sr/(1-sr^2)^(3/2)*(acos(sr)-sr*(1-sr^2)^(1/2));
elseif sr==1 % sphere
    error('use fun_P_sphere_iso instead')
end

% computation of elastic constants of matrix material
D0=inv(C0);
E0=1/D0(1,1);
mu0=0.5*C0(4,4);
nu0=E0/(2*mu0)-1;

% Eshelby Tensor components
S1111=3/(8*(1-nu0))*sr^2/(sr^2-1)+1/(4*(1-nu0))*(1-2*nu0-9/(4*(sr^2-1)))*g;
S2222=S1111;
S3333=1/(2*(1-nu0))*(1-2*nu0+(3*sr^2-1)/(sr^2-1)-(1-2*nu0+3*sr^2/(sr^2-1))*g);
S1122=1/(4*(1-nu0))*(sr^2/(2*(sr^2-1))-(1-2*nu0+3/(4*(sr^2-1)))*g);
S2211=S1122;
S1133=-1/(2*(1-nu0))*sr^2/(sr^2-1)+1/(4*(1-nu0))*(3*sr^2/(sr^2-1)-(1-2*nu0))*g;
S2233=S1133;
S3311=1/(2*(1-nu0))*(-1+2*nu0+1/(1-sr^2))+1/(4*(1-nu0))*(2*(1-2*nu0)-3/(1-sr^2))*g;
S3322=S3311;
S1212=1/(8*(1-nu0))*sr^2/(sr^2-1)+1/(4*(1-nu0))*(1-2*nu0-3/(4*(sr^2-1)))*g;
S3131=1/(4*(1-nu0))*(1-2*nu0-(sr^2+1)/(sr^2-1))-1/(8*(1-nu0))*(1-2*nu0-3*(sr^2+1)/(sr^2-1))*g;
S2323=S3131;

S  = [1*S1111,1*S1122,1*S1133,0,0,0;...
    1*S2211,1*S2222,1*S2233,0,0,0;...
    1*S3311,1*S3322,1*S3333,0,0,0;...
    0,0,0,2*S2323,0,0;...
    0,0,0,0,2*S3131,0;...
    0,0,0,0,0,2*S1212];

% Hill Tensor
P=S/C0;
