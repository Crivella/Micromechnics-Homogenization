% function for computation of P-tensor for ellipsoidal inclusions in 
% transversely isotropic matrix material with axis 3 as axis of rotational 
% symmetry
%
%-----------------------------------------------------------------------
% INPUT: 
% C0       ... (trans. iso.) matrix stiffness tensor C0
% asp   ... aspect ratio (a1/a2)
% slend ... slenderness ratio (a1/a3)
% toli     ... absolute precision, default=1e-9
%
% OUTPUT: P-tensor
% Pellips ... P-tensor
%
%-----------------------------------------------------------------------
%
function[Pellips]=fun_P_spheroid_transiso(C0,asp,slend,toli)
% Precision
if nargin<3.5
    toli=6;
end
mintol=-1;
maxtol=-16;
tol_list=logspace(mintol,maxtol,16);
abstol=tol_list(toli);
reltol=abstol;

% determination of (independent) components of stiffness tensor C0
C1111=C0(1,1);
C3333=C0(3,3);
C1122=C0(1,2);
C2233=C0(2,3);
C2323=C0(4,4)/2;

funP1111 = @(phi,theta) (4..*cos(phi).^2.*sin(theta).^3.*(C2323.*C3333.*slend.^4.*cos(theta).^4 + slend.^2.*cos(theta).^2.*((C2323.^2 + 0.5.*C1111.*C3333 - 0.5.*C1122.*C3333).* cos(phi).^2 + asp.^2.*(-1..*C2233.^2 - 2..*C2233.*C2323 + C1111.*C3333).* sin(phi).^2).*sin(theta).^2 + C2323.*((0.5.*C1111 - 0.5.*C1122).* cos(phi).^4 + asp.^2.*(1.5.*C1111 - 0.5.*C1122).*cos(phi).^2.*sin(phi).^2 + asp.^4.*C1111.*sin(phi).^4).*sin(theta).^4))./ (C2323.^2.*C3333.*slend.^6.*cos(theta).^6 - 2..*C2323.*(0.5.*C2233.^2 + 1..*C2233.*C2323 - 0.75.*C1111.*C3333 + 0.25.*C1122.*C3333).*slend.^4.*cos(theta).^4.*(cos(phi).^2 + asp.^2.*sin(phi).^2).*sin(theta).^2 + slend.^2.*cos(theta).^2.* ((C1122.*C2233.*(0.5.*C2233 + 1..*C2323) + 0.5.*C1111.^2.*C3333 + C1111.*(-0.5.*C2233.^2 - 1..*C2233.*C2323 + 1..*C2323.^2 - 0.5.*C1122.*C3333)).*cos(phi).^4 + 1..*asp.^2.* (C1122.*C2233.*(1..*C2233 + 2..*C2323) + 1..*C1111.^2.*C3333 + C1111.*(-1..*C2233.^2 - 2..*C2233.*C2323 + 2..*C2323.^2 - 1..*C1122.*C3333)).* cos(phi).^2.*sin(phi).^2 + asp.^4.*(C1122.*C2233.*(0.5.*C2233 + 1..*C2323) + 0.5.*C1111.^2.*C3333 + C1111.*(-0.5.*C2233.^2 - 1..*C2233.*C2323 + 1..*C2323.^2 - 0.5.*C1122.*C3333)).*sin(phi).^4).*sin(theta).^4 + 0.5.*C1111.*(1..*C1111 - 1..*C1122).*C2323.*(cos(phi).^6 + 3..*asp.^2.*cos(phi).^4.*sin(phi).^2 + 3..*asp.^4.*cos(phi).^2.*sin(phi).^4 + asp.^6.*sin(phi).^6).*sin(theta).^6);
funP1122 = @(phi,theta) (4..*asp.^2.*cos(phi).^2.*sin(phi).^2.*sin(theta).^5.* ((C2233.^2 + 2..*C2233.*C2323 + C2323.^2 - 0.5.*C1111.*C3333 - 0.5.*C1122.*C3333).*slend.^2.*cos(theta).^2 + (-0.5.*C1111 - 0.5.*C1122).* C2323.*(cos(phi).^2 + asp.^2.*sin(phi).^2).*sin(theta).^2))./ (C2323.^2.*C3333.*slend.^6.*cos(theta).^6 - 2..*C2323.*(0.5.*C2233.^2 + 1..*C2233.*C2323 - 0.75.*C1111.*C3333 + 0.25.*C1122.*C3333).*slend.^4.*cos(theta).^4.*(cos(phi).^2 + asp.^2.*sin(phi).^2).*sin(theta).^2 + slend.^2.*cos(theta).^2.* ((C1122.*C2233.*(0.5.*C2233 + 1..*C2323) + 0.5.*C1111.^2.*C3333 + C1111.*(-0.5.*C2233.^2 - 1..*C2233.*C2323 + 1..*C2323.^2 - 0.5.*C1122.*C3333)).*cos(phi).^4 + 1..*asp.^2.* (C1122.*C2233.*(1..*C2233 + 2..*C2323) + 1..*C1111.^2.*C3333 + C1111.*(-1..*C2233.^2 - 2..*C2233.*C2323 + 2..*C2323.^2 - 1..*C1122.*C3333)).* cos(phi).^2.*sin(phi).^2 + asp.^4.*(C1122.*C2233.*(0.5.*C2233 + 1..*C2323) + 0.5.*C1111.^2.*C3333 + C1111.*(-0.5.*C2233.^2 - 1..*C2233.*C2323 + 1..*C2323.^2 - 0.5.*C1122.*C3333)).*sin(phi).^4).*sin(theta).^4 + 0.5.*C1111.*(1..*C1111 - 1..*C1122).*C2323.*(cos(phi).^6 + 3..*asp.^2.*cos(phi).^4.*sin(phi).^2 + 3..*asp.^4.*cos(phi).^2.*sin(phi).^4 + asp.^6.*sin(phi).^6).*sin(theta).^6);
funP1133 = @(phi,theta) (4..*slend.^2.*cos(phi).^2.*cos(theta).^2.*sin(theta).^3.* ((-1..*C2233 - 1..*C2323).*C2323.*slend.^2.*cos(theta).^2 + (-0.5.*C1111.*C2233 + 0.5.*C1122.*C2233 - 0.5.*C1111.*C2323 + 0.5.*C1122.*C2323).*(cos(phi).^2 + asp.^2.*sin(phi).^2).*sin(theta).^2))./ (C2323.^2.*C3333.*slend.^6.*cos(theta).^6 - 2..*C2323.*(0.5.*C2233.^2 + 1..*C2233.*C2323 - 0.75.*C1111.*C3333 + 0.25.*C1122.*C3333).*slend.^4.*cos(theta).^4.*(cos(phi).^2 + asp.^2.*sin(phi).^2).*sin(theta).^2 + slend.^2.*cos(theta).^2.* ((C1122.*C2233.*(0.5.*C2233 + 1..*C2323) + 0.5.*C1111.^2.*C3333 + C1111.*(-0.5.*C2233.^2 - 1..*C2233.*C2323 + 1..*C2323.^2 - 0.5.*C1122.*C3333)).*cos(phi).^4 + 1..*asp.^2.* (C1122.*C2233.*(1..*C2233 + 2..*C2323) + 1..*C1111.^2.*C3333 + C1111.*(-1..*C2233.^2 - 2..*C2233.*C2323 + 2..*C2323.^2 - 1..*C1122.*C3333)).* cos(phi).^2.*sin(phi).^2 + asp.^4.*(C1122.*C2233.*(0.5.*C2233 + 1..*C2323) + 0.5.*C1111.^2.*C3333 + C1111.*(-0.5.*C2233.^2 - 1..*C2233.*C2323 + 1..*C2323.^2 - 0.5.*C1122.*C3333)).*sin(phi).^4).*sin(theta).^4 + 0.5.*C1111.*(1..*C1111 - 1..*C1122).*C2323.*(cos(phi).^6 + 3..*asp.^2.*cos(phi).^4.*sin(phi).^2 + 3..*asp.^4.*cos(phi).^2.*sin(phi).^4 + asp.^6.*sin(phi).^6).*sin(theta).^6);
funP3333 = @(phi,theta) (4..*slend.^2.*cos(theta).^2.*sin(theta).*(C2323.^2.*slend.^4.*cos(theta).^4 + 1.5.*(1..*C1111 - 1./3.*C1122).*C2323.*slend.^2.*cos(theta).^2.* (cos(phi).^2 + asp.^2.*sin(phi).^2).*sin(theta).^2 + 0.5.*C1111.*(1..*C1111 - 1..*C1122).*(cos(phi).^4 + 2..*asp.^2.*cos(phi).^2.* sin(phi).^2 + asp.^4.*sin(phi).^4).*sin(theta).^4))./ (C2323.^2.*C3333.*slend.^6.*cos(theta).^6 - 2..*C2323.*(0.5.*C2233.^2 + 1..*C2233.*C2323 - 0.75.*C1111.*C3333 + 0.25.*C1122.*C3333).*slend.^4.*cos(theta).^4.*(cos(phi).^2 + asp.^2.*sin(phi).^2).*sin(theta).^2 + slend.^2.*cos(theta).^2.* ((C1122.*C2233.*(0.5.*C2233 + 1..*C2323) + 0.5.*C1111.^2.*C3333 + C1111.*(-0.5.*C2233.^2 - 1..*C2233.*C2323 + 1..*C2323.^2 - 0.5.*C1122.*C3333)).*cos(phi).^4 + 1..*asp.^2.* (C1122.*C2233.*(1..*C2233 + 2..*C2323) + 1..*C1111.^2.*C3333 + C1111.*(-1..*C2233.^2 - 2..*C2233.*C2323 + 2..*C2323.^2 - 1..*C1122.*C3333)).* cos(phi).^2.*sin(phi).^2 + asp.^4.*(C1122.*C2233.*(0.5.*C2233 + 1..*C2323) + 0.5.*C1111.^2.*C3333 + C1111.*(-0.5.*C2233.^2 - 1..*C2233.*C2323 + 1..*C2323.^2 - 0.5.*C1122.*C3333)).*sin(phi).^4).*sin(theta).^4 + 0.5.*C1111.*(1..*C1111 - 1..*C1122).*C2323.*(cos(phi).^6 + 3..*asp.^2.*cos(phi).^4.*sin(phi).^2 + 3..*asp.^4.*cos(phi).^2.*sin(phi).^4 + asp.^6.*sin(phi).^6).*sin(theta).^6);
funP2323 = @(phi,theta) (C2323.*C3333.*slend.^6.*cos(theta).^6.*sin(theta) - 2..*slend.^4.*cos(theta).^4.*((0.5.*C2233.^2 + 1..*C2233.*C2323 - 0.5.*C1111.*C3333).*cos(phi).^2 + asp.^2.*(1..*C2233.*C2323 - 0.25.*C1111.*C3333 + 0.25.*C1122.*C3333).*sin(phi).^2).*sin(theta).^3 + 1..*slend.^2.*cos(theta).^2.*(1..*C1111.*C2323.*cos(phi).^4 - 1..*asp.^2.*(1..*C1111.*C2233 - 1..*C1122.*C2233 - 2..*C1111.*C2323).* cos(phi).^2.*sin(phi).^2 + asp.^4.*(-1..*C1111.*C2233 + 1..*C1122.*C2233 + 1..*C1111.*C2323).*sin(phi).^4).*sin(theta).^5 + 0.5.*asp.^2.*C1111.*(1..*C1111 - 1..*C1122).*sin(phi).^2.* (cos(phi).^4 + 2..*asp.^2.*cos(phi).^2.*sin(phi).^2 + asp.^4.*sin(phi).^4).* sin(theta).^7)./(C2323.^2.*C3333.*slend.^6.*cos(theta).^6 - 2..*C2323.*(0.5.*C2233.^2 + 1..*C2233.*C2323 - 0.75.*C1111.*C3333 + 0.25.*C1122.*C3333).*slend.^4.*cos(theta).^4.*(cos(phi).^2 + asp.^2.*sin(phi).^2).*sin(theta).^2 + slend.^2.*cos(theta).^2.* ((C1122.*C2233.*(0.5.*C2233 + 1..*C2323) + 0.5.*C1111.^2.*C3333 + C1111.*(-0.5.*C2233.^2 - 1..*C2233.*C2323 + 1..*C2323.^2 - 0.5.*C1122.*C3333)).*cos(phi).^4 + 1..*asp.^2.* (C1122.*C2233.*(1..*C2233 + 2..*C2323) + 1..*C1111.^2.*C3333 + C1111.*(-1..*C2233.^2 - 2..*C2233.*C2323 + 2..*C2323.^2 - 1..*C1122.*C3333)).* cos(phi).^2.*sin(phi).^2 + asp.^4.*(C1122.*C2233.*(0.5.*C2233 + 1..*C2323) + 0.5.*C1111.^2.*C3333 + C1111.*(-0.5.*C2233.^2 - 1..*C2233.*C2323 + 1..*C2323.^2 - 0.5.*C1122.*C3333)).*sin(phi).^4).*sin(theta).^4 + 0.5.*C1111.*(1..*C1111 - 1..*C1122).*C2323.*(cos(phi).^6 + 3..*asp.^2.*cos(phi).^4.*sin(phi).^2 + 3..*asp.^4.*cos(phi).^2.*sin(phi).^4 + asp.^6.*sin(phi).^6).*sin(theta).^6);
funP1212 = @(phi,theta) (sin(theta).^3.*(C2323.*C3333.*slend.^4.*cos(theta).^4.* (cos(phi).^2 + asp.^2.*sin(phi).^2) + 4..*slend.^2.*cos(theta).^2.* ((-0.25.*C2233.^2 - 0.5.*C2233.*C2323 + 0.25.*C1111.*C3333).*cos(phi).^4 + 0.5.*asp.^2.*(1..*C2233.^2 + 2..*C2233.*C2323 + 2..*C2323.^2 - 1..*C1122.*C3333).*cos(phi).^2.*sin(phi).^2 + asp.^4.*(-0.25.*C2233.^2 - 0.5.*C2233.*C2323 + 0.25.*C1111.*C3333).* sin(phi).^4).*sin(theta).^2 + 1..*C2323.*(1..*C1111.*cos(phi).^6 + 1..*asp.^2.*(1..*C1111 - 2..*C1122).*cos(phi).^4.*sin(phi).^2 + 1..*asp.^4.*(1..*C1111 - 2..*C1122).*cos(phi).^2.*sin(phi).^4 + 1..*asp.^6.*C1111.*sin(phi).^6).*sin(theta).^4))./ (C2323.^2.*C3333.*slend.^6.*cos(theta).^6 - 2..*C2323.*(0.5.*C2233.^2 + 1..*C2233.*C2323 - 0.75.*C1111.*C3333 + 0.25.*C1122.*C3333).*slend.^4.*cos(theta).^4.*(cos(phi).^2 + asp.^2.*sin(phi).^2).*sin(theta).^2 + slend.^2.*cos(theta).^2.* ((C1122.*C2233.*(0.5.*C2233 + 1..*C2323) + 0.5.*C1111.^2.*C3333 + C1111.*(-0.5.*C2233.^2 - 1..*C2233.*C2323 + 1..*C2323.^2 - 0.5.*C1122.*C3333)).*cos(phi).^4 + 1..*asp.^2.* (C1122.*C2233.*(1..*C2233 + 2..*C2323) + 1..*C1111.^2.*C3333 + C1111.*(-1..*C2233.^2 - 2..*C2233.*C2323 + 2..*C2323.^2 - 1..*C1122.*C3333)).* cos(phi).^2.*sin(phi).^2 + asp.^4.*(C1122.*C2233.*(0.5.*C2233 + 1..*C2323) + 0.5.*C1111.^2.*C3333 + C1111.*(-0.5.*C2233.^2 - 1..*C2233.*C2323 + 1..*C2323.^2 - 0.5.*C1122.*C3333)).*sin(phi).^4).*sin(theta).^4 + 0.5.*C1111.*(1..*C1111 - 1..*C1122).*C2323.*(cos(phi).^6 + 3..*asp.^2.*cos(phi).^4.*sin(phi).^2 + 3..*asp.^4.*cos(phi).^2.*sin(phi).^4 + asp.^6.*sin(phi).^6).*sin(theta).^6);

% computation of tensor components by integration over unit sphere
P1111=4/(16*pi)*integral2(funP1111,0,pi/2,0,pi,'abstol',abstol,'reltol',reltol);
P1122=4/(16*pi)*integral2(funP1122,0,pi/2,0,pi,'abstol',abstol,'reltol',reltol);
P1133=4/(16*pi)*integral2(funP1133,0,pi/2,0,pi,'abstol',abstol,'reltol',reltol);
P2222=P1111;
P2233=P1133;
P3333=4/(16*pi)*integral2(funP3333,0,pi/2,0,pi,'abstol',abstol,'reltol',reltol);
P2323=4/(16*pi)*integral2(funP2323,0,pi/2,0,pi,'abstol',abstol,'reltol',reltol);
P1313=P2323;
P1212=4/(16*pi)*integral2(funP1212,0,pi/2,0,pi,'abstol',abstol,'reltol',reltol);

% compilation of P-tensor in compressed notation
Pellips=zeros(6,6);
Pellips(1,1)=P1111;
Pellips(2,2)=P2222;
Pellips(3,3)=P3333;
Pellips(1,2)=P1122;
Pellips(2,1)=Pellips(1,2);
Pellips(1,3)=P1133;
Pellips(3,1)=Pellips(1,3);
Pellips(2,3)=P2233;
Pellips(3,2)=Pellips(2,3);

Pellips(4,4)=2*P2323;
Pellips(5,5)=2*P1313;
Pellips(6,6)=2*P1212;
