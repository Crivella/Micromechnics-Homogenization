function [comp] = prep_comp_AI(comp)
% %load
load('fibprops.mat','rho');

% consider measured crystallinity
CI=comp.fib.CI./(comp.fib.CI+(1-comp.fib.CI)*rho.amcel/rho.crycel);

% consider measured chemical composition
fld=fieldnames(comp.fib.mchem);
tmpsum1=0;
for i=1:numel(fld)
    if ~strcmp(fld{i},'extr')
        tmpsum1=tmpsum1+comp.fib.mchem.(fld{i});
    end
end
comp.fib.mchem.extr=1-tmpsum1;
tmpsum1=comp.fib.mchem.extr+tmpsum1;
for i=1:numel(fld)
    comp.fib.mchem.(fld{i})=comp.fib.mchem.(fld{i})/tmpsum1;
end

% crycel+amcel
comp.fib.mchem.crycel=CI*comp.fib.mchem.totcel;
comp.fib.mchem.amcel=(1-CI)*comp.fib.mchem.totcel;
comp.fib.mchem=rmfield(comp.fib.mchem,'totcel');

% volumes
fld=fieldnames(comp.fib.mchem);
for i=1:numel(fld)
    comp.fib.vchem.(fld{i})=comp.fib.mchem.(fld{i})/rho.(fld{i});
end

% volume fractions
tmpsum=0;
for i=1:numel(fld)
    tmpsum=tmpsum+comp.fib.vchem.(fld{i});
end
tmpsum1=0;
for i=1:numel(fld)
    f.(fld{i})=comp.fib.vchem.(fld{i})./tmpsum;
    tmpsum1=tmpsum1+f.(fld{i});
end
%end

%B)composite-related fractions
if isnan(comp.mat.fair)
    v.air=0;
else
    v.air=comp.mat.fair;
end
v.fib=(1-v.air)/(1+(1-comp.fib.vfrac)/comp.fib.vfrac);
v.lum=comp.fib.lumpor*v.fib;
v.cw=v.fib*(1-comp.fib.lumpor);
v.mat=1-v.fib-v.air;
fld=fieldnames(f);
for i=1:numel(fld)
    v.(fld{i})=f.(fld{i})*v.cw;
end
comp.v=v;
end

