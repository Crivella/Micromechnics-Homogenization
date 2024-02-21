function [comp] = prep_comp_noexp(comp)
%load
load('fibprops.mat') %load fiber database
load('matprops.mat') % load matrix database

% set chemophysical fiber properties, if not defined
physid=find(strcmp(comp.fib.name,cellstr(phys.name)));
chemid=find(strcmp(comp.fib.name,cellstr(chem.name)));
%stdfib.rho=phys.rho(physid,2);
%stdfib.lumpor=phys.lumpor(physid,2);
comp.fib.rhocalc=phys.rhocalc(physid);
fld={'l','d','rho','MFA'};
for i=1:numel(fld)
    if isnan(comp.fib.(fld{i}))
        comp.fib.(fld{i})=phys.(fld{i})(physid,comp.case);
    end
end
if isnan(comp.fib.lumpor)
    comp.fib.lumpor=max(0,1-comp.fib.rho./phys.rhocalc(physid));
end

if isnan(comp.fib.ar) % if NaN, calculate based on l and d
    comp.fib.ar=comp.fib.l/comp.fib.d*1000;
end
if isnan(comp.fib.ar) % if still NaN, take standard ar
    comp.fib.ar=phys.ar(physid);
end

% set matrix properties, if not defined
fld={'rho','E','nu'};
matdid=find(strcmp(comp.mat.name,stdmat.name));
for i=1:numel(fld)
    if isnan(comp.mat.(fld{i}))
        comp.mat.(fld{i})=stdmat.(fld{i})(matdid);
    end
end

% % set experimental properties to zero, if not defined
% fld=fieldnames(comp.exp);
% for i=1:numel(fld)
%     if isnan(comp.exp.(fld{i}))
%         comp.exp.(fld{i})=0;
%     end
% end

% A)fiber-related fractions
if isnan(comp.fib.mchem.totcel) && isnan(comp.fib.CI)
    % CASE NO CHEMISTRY MEASURED
    fld=fieldnames(chemv);
    for i=1:numel(fld)
        f.(fld{i})=chemv.(fld{i})(chemid,comp.case);
    end
else
    if ~isnan(comp.fib.CI)
        %consider measured crystallinity
        CI=comp.fib.CI./(comp.fib.CI+(1-comp.fib.CI)*rho.amcel/rho.crycel);
    else
        %consider standard crystallinity
        CI=chem.crysm(chemid,comp.case);
    end
    if ~isnan(comp.fib.mchem.totcel)
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
    else
        % consider standard chemical composition
        switch comp.case
            case 1; tmpchem=chemmin.mn;
            case 2; tmpchem=chemmax.mn;
            case 3; tmpchem=chem.mn;
        end
        comp.fib.mchem.totcel=(tmpchem.crycel(chemid)+tmpchem.amcel(chemid))/100;
        fld=fieldnames(chem.m);
        for i=1:numel(fld)
            if ~strcmp(fld{i},'totcel')
                comp.fib.mchem.(fld{i})=tmpchem.(fld{i})(chemid)/100;
            end
        end
    end
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
end

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

