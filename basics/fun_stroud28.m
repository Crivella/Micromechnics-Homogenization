% STROUT POINT, 28 Point integration
% last modified: 21.9.2017
% to do: ---
function [stroud] = fun_stroud28()
% Stroud Points for 56 point integration
r28=sqrt((9 -4*sqrt(3))/33);
s28=sqrt((15+8*sqrt(3))/33);
t28=sqrt(1/3);
u28=sqrt((15-8*sqrt(3))/33);
v28=sqrt((9 +4*sqrt(3))/33);

z28=[t28 t28 t28 t28 r28 r28 r28 r28 r28 r28 r28 r28 s28 s28 v28 v28 v28 v28 v28 u28 u28 u28 u28 s28 s28 v28 v28 v28];
y28=[t28 -t28 t28 -t28 r28 -r28 r28 -r28 s28 -s28 s28 -s28 r28 -r28 -v28 u28 -u28 u28 -u28 v28 -v28 v28 -v28 r28 -r28 v28 -v28 v28];
x28=[t28 t28 -t28 -t28 s28 s28 -s28 -s28 r28 r28 -r28 -r28 r28 r28 -u28 v28 v28 -v28 -v28 v28 v28 -v28 -v28 -r28 -r28 u28 u28 -u28];

a28=9/280;
b28=(122+9*sqrt(3))/3360;
c28=(122-9*sqrt(3))/3360;
stroud.weight=[a28 a28 a28 a28 b28 b28 b28 b28 b28 b28 b28 b28 b28 b28 c28 c28 c28 c28 c28 c28 c28 c28 c28 b28 b28 c28 c28 c28];

[stroud.azi,elev,~]=cart2sph(x28,y28,z28);
stroud.zeni=pi/2-elev;
end
