% function for computation of R-tensor related to "slightly" imperfect
% interfaces for ellisoidal inclusions in matrix-inclusion problems, 
% see Qu 1993a,b and Dinzart & Sabar 2017
%
%-----------------------------------------------------------------------
% INPUT: 
% asp   . aspect ratio (a1./a2)
% slend . slenderness ratio (a1./a3)
% alpha . interface sliding parameter, alpha>=0 [m./GPa]
% beta . interface separation parameter, beta>=0 [m./GPa]
% tol . tolerance integer

%
% OUTPUT: R-tensor
% R .. R-tensor
%
%-----------------------------------------------------------------------
%
function[R]=fun_R_ellipsoidDinzart(asp,slend,alpha,beta,toli)
if nargin<4.5
    toli=6;
end
%toli=6;asp=1;slend=2;
mintol=-1;
maxtol=-16;
tol_list=logspace(mintol,maxtol,16);
abstol=tol_list(toli);
reltol=abstol;

%eta=[beta 0 0; 0 alpha 0; 0 0 alpha]; % interface stiffness in local (local ellipsoidal coord)
%R=zeros(3,3,3,3);
a1=1; % everything normalized
a2=a1/asp;
a3=a1/slend;
%Vi=a1.*a2.*a3;
%%
Pint1111= @(azi,zeni) 4.*sin(azi).^2.*cos(zeni).^2./a1.^2./(sin(azi).^2.*cos(zeni).^2./a1.^2+sin(azi).^2.*sin(zeni).^2./a2.^2+cos(azi).^2./a3.^2).^(1./2).*sin(zeni);
Pint1122 = 0;
Pint1133 = 0;
Pint1123 = 0;
Pint1131 = @(azi,zeni) 2.*sin(azi).*cos(zeni)./a1.*cos(azi)./a3./(sin(azi).^2.*cos(zeni).^2./a1.^2+sin(azi).^2.*sin(zeni).^2./a2.^2+cos(azi).^2./a3.^2).^(1./2).*sin(zeni); 
Pint1112 = @(azi,zeni) 2.*sin(azi).^2.*cos(zeni)./a1.*sin(zeni)./a2./(sin(azi).^2.*cos(zeni).^2./a1.^2+sin(azi).^2.*sin(zeni).^2./a2.^2+cos(azi).^2./a3.^2).^(1./2).*sin(zeni);
Pint2222 = @(azi,zeni) 4.*sin(azi).^2.*sin(zeni).^2./a2.^2./(sin(azi).^2.*cos(zeni).^2./a1.^2+sin(azi).^2.*sin(zeni).^2./a2.^2+cos(azi).^2./a3.^2).^(1./2).*sin(zeni);
Pint2233 = 0; 
Pint2223 = @(azi,zeni) 2.*sin(azi).*sin(zeni)./a2.*cos(azi)./a3./(sin(azi).^2.*cos(zeni).^2./a1.^2+sin(azi).^2.*sin(zeni).^2./a2.^2+cos(azi).^2./a3.^2).^(1./2).*sin(zeni);
Pint2231 = 0;
Pint2212 = @(azi,zeni) 2.*sin(azi).^2.*cos(zeni)./a1.*sin(zeni)./a2./(sin(azi).^2.*cos(zeni).^2./a1.^2+sin(azi).^2.*sin(zeni).^2./a2.^2+cos(azi).^2./a3.^2).^(1./2).*sin(zeni);
Pint3333 = @(azi,zeni) 4.*cos(azi).^2./a3.^2./(sin(azi).^2.*cos(zeni).^2./a1.^2+sin(azi).^2.*sin(zeni).^2./a2.^2+cos(azi).^2./a3.^2).^(1./2).*sin(zeni);
Pint3323 = @(azi,zeni) 2.*sin(azi).*sin(zeni)./a2.*cos(azi)./a3./(sin(azi).^2.*cos(zeni).^2./a1.^2+sin(azi).^2.*sin(zeni).^2./a2.^2+cos(azi).^2./a3.^2).^(1./2).*sin(zeni);
Pint3331 = @(azi,zeni) 2.*sin(azi).*cos(zeni)./a1.*cos(azi)./a3./(sin(azi).^2.*cos(zeni).^2./a1.^2+sin(azi).^2.*sin(zeni).^2./a2.^2+cos(azi).^2./a3.^2).^(1./2).*sin(zeni);
Pint3312 = 0;
Pint2323 = @(azi,zeni) (cos(azi).^2./a3.^2+sin(azi).^2.*sin(zeni).^2./a2.^2)./(sin(azi).^2.*cos(zeni).^2./a1.^2+sin(azi).^2.*sin(zeni).^2./a2.^2+cos(azi).^2./a3.^2).^(1./2).*sin(zeni);
Pint2331 = @(azi,zeni) sin(azi).^2.*cos(zeni)./a1.*sin(zeni)./a2./(sin(azi).^2.*cos(zeni).^2./a1.^2+sin(azi).^2.*sin(zeni).^2./a2.^2+cos(azi).^2./a3.^2).^(1./2).*sin(zeni);
Pint2312 = @(azi,zeni) sin(azi).*cos(zeni)./a1.*cos(azi)./a3./(sin(azi).^2.*cos(zeni).^2./a1.^2+sin(azi).^2.*sin(zeni).^2./a2.^2+cos(azi).^2./a3.^2).^(1./2).*sin(zeni);
Pint3131 = @(azi,zeni) (cos(azi).^2./a3.^2+sin(azi).^2.*cos(zeni).^2./a1.^2)./(sin(azi).^2.*cos(zeni).^2./a1.^2+sin(azi).^2.*sin(zeni).^2./a2.^2+cos(azi).^2./a3.^2).^(1./2).*sin(zeni);
Pint3112 = @(azi,zeni) sin(azi).*sin(zeni)./a2.*cos(azi)./a3./(sin(azi).^2.*cos(zeni).^2./a1.^2+sin(azi).^2.*sin(zeni).^2./a2.^2+cos(azi).^2./a3.^2).^(1./2).*sin(zeni); 
Pint1212 = @(azi,zeni) (sin(azi).^2.*sin(zeni).^2./a2.^2+sin(azi).^2.*cos(zeni).^2./a1.^2)./(sin(azi).^2.*cos(zeni).^2./a1.^2+sin(azi).^2.*sin(zeni).^2./a2.^2+cos(azi).^2./a3.^2).^(1./2).*sin(zeni);


Qint1111 = @(azi,zeni) sin(azi).^4.*cos(zeni).^4./a1.^4./(sin(azi).^2.*cos(zeni).^2./a1.^2+sin(azi).^2.*sin(zeni).^2./a2.^2+cos(azi).^2./a3.^2).^(3./2).*sin(zeni);
Qint1122 = @(azi,zeni) sin(azi).^4.*cos(zeni).^2./a1.^2.*sin(zeni).^2./a2.^2./(sin(azi).^2.*cos(zeni).^2./a1.^2+sin(azi).^2.*sin(zeni).^2./a2.^2+cos(azi).^2./a3.^2).^(3./2).*sin(zeni);
Qint1133 = @(azi,zeni) sin(azi).^2.*cos(zeni).^2./a1.^2.*cos(azi).^2./a3.^2./(sin(azi).^2.*cos(zeni).^2./a1.^2+sin(azi).^2.*sin(zeni).^2./a2.^2+cos(azi).^2./a3.^2).^(3./2).*sin(zeni);
Qint1123 = @(azi,zeni) sin(azi).^3.*cos(zeni).^2./a1.^2.*sin(zeni)./a2.*cos(azi)./a3./(sin(azi).^2.*cos(zeni).^2./a1.^2+sin(azi).^2.*sin(zeni).^2./a2.^2+cos(azi).^2./a3.^2).^(3./2).*sin(zeni);
Qint1131 = @(azi,zeni) sin(azi).^3.*cos(zeni).^3./a1.^3.*cos(azi)./a3./(sin(azi).^2.*cos(zeni).^2./a1.^2+sin(azi).^2.*sin(zeni).^2./a2.^2+cos(azi).^2./a3.^2).^(3./2).*sin(zeni);
Qint1112 = @(azi,zeni) sin(azi).^4.*cos(zeni).^3./a1.^3.*sin(zeni)./a2./(sin(azi).^2.*cos(zeni).^2./a1.^2+sin(azi).^2.*sin(zeni).^2./a2.^2+cos(azi).^2./a3.^2).^(3./2).*sin(zeni);
Qint2222 = @(azi,zeni) sin(azi).^4.*sin(zeni).^4./a2.^4./(sin(azi).^2.*cos(zeni).^2./a1.^2+sin(azi).^2.*sin(zeni).^2./a2.^2+cos(azi).^2./a3.^2).^(3./2).*sin(zeni);
Qint2233 = @(azi,zeni) sin(azi).^2.*sin(zeni).^2./a2.^2.*cos(azi).^2./a3.^2./(sin(azi).^2.*cos(zeni).^2./a1.^2+sin(azi).^2.*sin(zeni).^2./a2.^2+cos(azi).^2./a3.^2).^(3./2).*sin(zeni);
Qint2223 = @(azi,zeni) sin(azi).^3.*sin(zeni).^3./a2.^3.*cos(azi)./a3./(sin(azi).^2.*cos(zeni).^2./a1.^2+sin(azi).^2.*sin(zeni).^2./a2.^2+cos(azi).^2./a3.^2).^(3./2).*sin(zeni);
Qint2231 = @(azi,zeni) sin(azi).^3.*cos(zeni)./a1.*sin(zeni).^2./a2.^2.*cos(azi)./a3./(sin(azi).^2.*cos(zeni).^2./a1.^2+sin(azi).^2.*sin(zeni).^2./a2.^2+cos(azi).^2./a3.^2).^(3./2).*sin(zeni);
Qint2212 = @(azi,zeni) sin(azi).^4.*cos(zeni)./a1.*sin(zeni).^3./a2.^3./(sin(azi).^2.*cos(zeni).^2./a1.^2+sin(azi).^2.*sin(zeni).^2./a2.^2+cos(azi).^2./a3.^2).^(3./2).*sin(zeni);
Qint3333 = @(azi,zeni) cos(azi).^4./a3.^4./(sin(azi).^2.*cos(zeni).^2./a1.^2+sin(azi).^2.*sin(zeni).^2./a2.^2+cos(azi).^2./a3.^2).^(3./2).*sin(zeni);
Qint3323 = @(azi,zeni) sin(azi).*sin(zeni)./a2.*cos(azi).^3./a3.^3./(sin(azi).^2.*cos(zeni).^2./a1.^2+sin(azi).^2.*sin(zeni).^2./a2.^2+cos(azi).^2./a3.^2).^(3./2).*sin(zeni);
Qint3331 = @(azi,zeni) sin(azi).*cos(zeni)./a1.*cos(azi).^3./a3.^3./(sin(azi).^2.*cos(zeni).^2./a1.^2+sin(azi).^2.*sin(zeni).^2./a2.^2+cos(azi).^2./a3.^2).^(3./2).*sin(zeni);
Qint3312 = @(azi,zeni) sin(azi).^2.*cos(zeni)./a1.*sin(zeni)./a2.*cos(azi).^2./a3.^2./(sin(azi).^2.*cos(zeni).^2./a1.^2+sin(azi).^2.*sin(zeni).^2./a2.^2+cos(azi).^2./a3.^2).^(3./2).*sin(zeni);
Qint2323 = @(azi,zeni) sin(azi).^2.*sin(zeni).^2./a2.^2.*cos(azi).^2./a3.^2./(sin(azi).^2.*cos(zeni).^2./a1.^2+sin(azi).^2.*sin(zeni).^2./a2.^2+cos(azi).^2./a3.^2).^(3./2).*sin(zeni);
Qint2331 = @(azi,zeni) sin(azi).^2.*cos(zeni)./a1.*sin(zeni)./a2.*cos(azi).^2./a3.^2./(sin(azi).^2.*cos(zeni).^2./a1.^2+sin(azi).^2.*sin(zeni).^2./a2.^2+cos(azi).^2./a3.^2).^(3./2).*sin(zeni);
Qint2312 = @(azi,zeni) sin(azi).^3.*cos(zeni)./a1.*sin(zeni).^2./a2.^2.*cos(azi)./a3./(sin(azi).^2.*cos(zeni).^2./a1.^2+sin(azi).^2.*sin(zeni).^2./a2.^2+cos(azi).^2./a3.^2).^(3./2).*sin(zeni);
Qint3131 = @(azi,zeni) sin(azi).^2.*cos(zeni).^2./a1.^2.*cos(azi).^2./a3.^2./(sin(azi).^2.*cos(zeni).^2./a1.^2+sin(azi).^2.*sin(zeni).^2./a2.^2+cos(azi).^2./a3.^2).^(3./2).*sin(zeni);
Qint3112 = @(azi,zeni) sin(azi).^3.*cos(zeni).^2./a1.^2.*sin(zeni)./a2.*cos(azi)./a3./(sin(azi).^2.*cos(zeni).^2./a1.^2+sin(azi).^2.*sin(zeni).^2./a2.^2+cos(azi).^2./a3.^2).^(3./2).*sin(zeni);
Qint1212 = @(azi,zeni) sin(azi).^4.*cos(zeni).^2./a1.^2.*sin(zeni).^2./a2.^2./(sin(azi).^2.*cos(zeni).^2./a1.^2+sin(azi).^2.*sin(zeni).^2./a2.^2+cos(azi).^2./a3.^2).^(3./2).*sin(zeni);


R1111=alpha.*3./(16.*pi).*integral2(Pint1111,0,2.*pi,0,pi,'abstol',abstol,'reltol',reltol)+3./(4.*pi).*(beta-alpha).*integral2(Qint1111,0,2.*pi,0,pi,'abstol',abstol,'reltol',reltol);
R1122=3./(4.*pi).*(beta-alpha).*integral2(Qint1122,0,2.*pi,0,pi,'abstol',abstol,'reltol',reltol);
R1133=3./(4.*pi).*(beta-alpha).*integral2(Qint1133,0,2.*pi,0,pi,'abstol',abstol,'reltol',reltol);
R1123=3./(4.*pi).*(beta-alpha).*integral2(Qint1123,0,2.*pi,0,pi,'abstol',abstol,'reltol',reltol);
R1131=alpha.*3./(16.*pi).*integral2(Pint1131,0,2.*pi,0,pi,'abstol',abstol,'reltol',reltol)+3./(4.*pi).*(beta-alpha).*integral2(Qint1131,0,2.*pi,0,pi,'abstol',abstol,'reltol',reltol);
R1112=alpha.*3./(16.*pi).*integral2(Pint1112,0,2.*pi,0,pi,'abstol',abstol,'reltol',reltol)+3./(4.*pi).*(beta-alpha).*integral2(Qint1112,0,2.*pi,0,pi,'abstol',abstol,'reltol',reltol);
R2222=alpha.*3./(16.*pi).*integral2(Pint2222,0,2.*pi,0,pi,'abstol',abstol,'reltol',reltol)+3./(4.*pi).*(beta-alpha).*integral2(Qint2222,0,2.*pi,0,pi,'abstol',abstol,'reltol',reltol);
R2233=3./(4.*pi).*(beta-alpha).*integral2(Qint2233,0,2.*pi,0,pi,'abstol',abstol,'reltol',reltol);
R2223=alpha.*3./(16.*pi).*integral2(Pint2223,0,2.*pi,0,pi,'abstol',abstol,'reltol',reltol)+3./(4.*pi).*(beta-alpha).*integral2(Qint2223,0,2.*pi,0,pi,'abstol',abstol,'reltol',reltol);
R2231=3./(4.*pi).*(beta-alpha).*integral2(Qint2231,0,2.*pi,0,pi,'abstol',abstol,'reltol',reltol);
R2212=alpha.*3./(16.*pi).*integral2(Pint2212,0,2.*pi,0,pi,'abstol',abstol,'reltol',reltol)+3./(4.*pi).*(beta-alpha).*integral2(Qint2212,0,2.*pi,0,pi,'abstol',abstol,'reltol',reltol);
R3333=alpha.*3./(16.*pi).*integral2(Pint3333,0,2.*pi,0,pi,'abstol',abstol,'reltol',reltol)+3./(4.*pi).*(beta-alpha).*integral2(Qint3333,0,2.*pi,0,pi,'abstol',abstol,'reltol',reltol);
R3323=alpha.*3./(16.*pi).*integral2(Pint3323,0,2.*pi,0,pi,'abstol',abstol,'reltol',reltol)+3./(4.*pi).*(beta-alpha).*integral2(Qint3323,0,2.*pi,0,pi,'abstol',abstol,'reltol',reltol);
R3331=alpha.*3./(16.*pi).*integral2(Pint3331,0,2.*pi,0,pi,'abstol',abstol,'reltol',reltol)+3./(4.*pi).*(beta-alpha).*integral2(Qint3331,0,2.*pi,0,pi,'abstol',abstol,'reltol',reltol);
R3312=3./(4.*pi).*(beta-alpha).*integral2(Qint3312,0,2.*pi,0,pi,'abstol',abstol,'reltol',reltol);
R2323=alpha.*3./(16.*pi).*integral2(Pint2323,0,2.*pi,0,pi,'abstol',abstol,'reltol',reltol)+3./(4.*pi).*(beta-alpha).*integral2(Qint2323,0,2.*pi,0,pi,'abstol',abstol,'reltol',reltol);
R2331=alpha.*3./(16.*pi).*integral2(Pint2331,0,2.*pi,0,pi,'abstol',abstol,'reltol',reltol)+3./(4.*pi).*(beta-alpha).*integral2(Qint2331,0,2.*pi,0,pi,'abstol',abstol,'reltol',reltol);
R2312=alpha.*3./(16.*pi).*integral2(Pint2312,0,2.*pi,0,pi,'abstol',abstol,'reltol',reltol)+3./(4.*pi).*(beta-alpha).*integral2(Qint2312,0,2.*pi,0,pi,'abstol',abstol,'reltol',reltol);
R3131=alpha.*3./(16.*pi).*integral2(Pint3131,0,2.*pi,0,pi,'abstol',abstol,'reltol',reltol)+3./(4.*pi).*(beta-alpha).*integral2(Qint3131,0,2.*pi,0,pi,'abstol',abstol,'reltol',reltol);
R3112=alpha.*3./(16.*pi).*integral2(Pint3112,0,2.*pi,0,pi,'abstol',abstol,'reltol',reltol)+3./(4.*pi).*(beta-alpha).*integral2(Qint3112,0,2.*pi,0,pi,'abstol',abstol,'reltol',reltol);
R1212=alpha.*3./(16.*pi).*integral2(Pint1212,0,2.*pi,0,pi,'abstol',abstol,'reltol',reltol)+3./(4.*pi).*(beta-alpha).*integral2(Qint1212,0,2.*pi,0,pi,'abstol',abstol,'reltol',reltol);

R=[R1111,R1122,R1133,sqrt(2).*R1123,sqrt(2).*R1131,sqrt(2).*R1112;
    0,R2222,R2233,sqrt(2).*R2223,sqrt(2).*R2231,sqrt(2).*R2212;
    0,0,R3333,sqrt(2).*R3323,sqrt(2).*R3331,sqrt(2).*R3312;
    0,0,0,2.*R2323,2.*R2331,2.*R2312;
    0,0,0,0,2.*R3131,2.*R3112;
    0,0,0,0,0,2.*R1212];
R=R+transpose(R)-R.*eye(6);

end
