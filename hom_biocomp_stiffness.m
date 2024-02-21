function [RVE,Phase] = hom_biocomp_stiffness(comp,tol)
if nargin<2
    tol=8;
end
%% elastic phase props
biocomp_materials

%% Homogenization
v=comp.v;
ar=1/comp.fib.ar(1); %cylinder=1e-20

% STEP1a: polymer network
Phase.hemcel=PhaseData('hemcel',Mat.hemcel,{'sphere'},v.hemcel);
Phase.lignin=PhaseData('lignin',Mat.lignin,{'sphere'},v.lignin);
Phase.pectin=PhaseData('pectin',Mat.pectin,{'sphere'},v.pectin);
Phase.extr=PhaseData('extr',Mat.pore,{'sphere'},v.extr+v.wax);
Phase.ash=PhaseData('ash',Mat.ash,{'sphere'},v.ash);

RVE.pn=RVEData1('pn',[Phase.hemcel,Phase.lignin,Phase.pectin,Phase.extr,Phase.ash],'SCS',{'matrix'},tol);
[RVE.pn,Phase.pn]=homE(RVE.pn);

% STEP1b: cellulose
Phase.amcel=PhaseData('amcel',Mat.amcel,{'matrix'},v.amcel);
Phase.crycel=PhaseData('crycel',Mat.crycel,{'spheroid',1e-20,[0,0]},v.crycel);

RVE.cel=RVEData1('cel',[Phase.amcel,Phase.crycel],'MT',{'spheroid',1e-20,[0,0]},tol);
[RVE.cel,Phase.cel]=homE(RVE.cel);

%% STEP2: cell wall
if comp.fib.MFA==0
    str=Phase.cel;
else
    nfam=20;
    str=repmat(Phase.cel,1,nfam);
    Mat.cel=MatData('cel','transiso',Phase.cel.mat.C,false,false);
    for i=1:nfam
        azi=(i-1)/nfam*2*pi;
        Phase.(['cel',num2str(i)])=PhaseData(['cel',num2str(i)],Mat.cel,{'spheroid',1e-20,[azi,comp.fib.MFA*pi/180]},Phase.cel.vol*1/nfam);
        str(i)=eval(['Phase.cel',num2str(i)]);
    end
end
%RVE.cw=RVEData1('fib',[str,Phase.pn],'MT',{'spheroid',comp.MLar,[0,0]},tol);
RVE.cw=RVEData1('cw',[str,Phase.pn],'MT',{'matrix'},tol);
[RVE.cw,Phase.cw]=homE(RVE.cw);
Phase.cw.mat.C=1/2*(Phase.cw.mat.C+transpose(Phase.cw.mat.C));% get something transversaly isotropic, get rid of strange assymetry !MORI_TANAKA PROBLEM!

%% STEP3: fiber bundle
Phase.lum=PhaseData('lum',Mat.pore,{'spheroid',comp.MLar,[0,0]},v.lum);
%Phase.ML=PhaseData('ml',Mat.lignin,{'matrix'},v.lignin);

RVE.fib=RVEData1('fib1',[Phase.lum,Phase.cw],'MT',{'spheroid',ar,'iso'},tol);
[RVE.fib,Phase.fib1]=homE(RVE.fib);

%% STEP4: biocomposite
Mat.mat=MatData('mat','iso',{'Enu',comp.mat.E,comp.mat.nu},false,false);
Phase.mat=PhaseData('mat',Mat.mat,{'matrix'},v.mat);
Phase.air=PhaseData('air',Mat.pore,{'sphere'},v.air);

% CONSIDER IF imperfection
%if ~isnan(comp.IF.mode) && comp.IF.mode~=0;    
    Phase.fib1.IF=comp.IF.par; 
%end

% CONSIDER ORIENTATION DISTRIBUTION
if strcmp(comp.fib.ori,'align') || strcmp(comp.fib.ori,'aligned') || strcmp(comp.fib.ori,'1D')
    Phase.fib1.shape.dis='single';
    if isnan(sum(comp.fib.oripar))
        Phase.fib1.shape.ori=[0,0];
    else
        Phase.fib1.shape.ori=[comp.fib.oripar(1)*pi/180,comp.fib.oripar(2)*pi/180];
    end
    str=Phase.fib1;
elseif strcmp(comp.fib.ori,'iso') || strcmp(comp.fib.ori,'3Diso')
    str=Phase.fib1;
elseif strcmp(comp.fib.ori,'2Diso') % in-plane 2D distribution
    nfam=20;
    str=repmat(Phase.fib1,1,nfam);
    Mat.fib=MatData('fib','transiso',Phase.fib1.mat.C,false,false);
    for i=1:nfam
        azi=(i-1)/nfam*pi;
        Phase.(['fib',num2str(i)])=PhaseData(['fib',num2str(i)],Mat.fib,{'spheroid',ar,[azi,pi/2]},v.fib*1/nfam);
        Phase.(['fib',num2str(i)]).IF=comp.IF.par;
        str(i)=eval(['Phase.fib',num2str(i)]);
    end
elseif strcmp(comp.fib.ori,'vMF') || strcmp(comp.fib.ori,'vMs')
    kappa=comp.fib.oripar(1);
    if kappa<2; nfam=14;
    elseif kappa<5; nfam=38;
    elseif kappa<10; nfam=86;
    elseif kappa<110; nfam=434;
    elseif kappa<1010; nfam=1202;
    elseif kappa<10010; nfam=5810;
        %if strcmp(comp.fib.ori,'vMF') && kappa<2; nfam=14;
        %elseif strcmp(comp.fib.ori,'vMF') && kappa<5; nfam=38;
        %elseif strcmp(comp.fib.ori,'vMF') && kappa<10; nfam=86;
        %elseif strcmp(comp.fib.ori,'vMF') && kappa<110; nfam=434;
        %elseif strcmp(comp.fib.ori,'vMF') && kappa<1010; nfam=1202;
        %elseif strcmp(comp.fib.ori,'vMF') && kappa<10010; nfam=5810;
        %elseif strcmp(comp.fib.ori,'vMs') && kappa<2; nfam=38;
        %elseif strcmp(comp.fib.ori,'vMs') && kappa<5; nfam=86;
        %elseif strcmp(comp.fib.ori,'vMs') && kappa<10; nfam=434;
        %elseif strcmp(comp.fib.ori,'vMs') && kappa<110; nfam=1202;
        %elseif strcmp(comp.fib.ori,'vMs') && kappa<1010; nfam=5810;
    else error('too many fibers required, reduce kappa')
    end
    str=repmat(Phase.fib1,1,nfam);
    Mat.fib=MatData('fib','transiso',Phase.fib1.mat.C,false,false);
    leb = getLebedevSphere(nfam);
    if strcmp(comp.fib.ori,'vMF');
        funact= kappa/sinh(kappa)/4/pi*exp(kappa*cos(leb.zeni));
    elseif strcmp(comp.fib.ori,'vMs')
        funvM = exp(kappa*cos(leb.zeni-pi/2))/(2*pi*besseli(0,kappa)); %vonMises
        leb1 = getLebedevSphere(5810);
        funvM1 = exp(kappa*cos(leb1.zeni-pi/2))/(2*pi*besseli(0,kappa)); %vonMises
        funact= funvM./sum(funvM1.*leb1.w); %vonMises normalized
    end
    funact(funact<0)=0;
    tmp=sum(funact.*leb.w);
    j=1;
    for i=1:nfam
        vol1=leb.w(i)*funact(i)/tmp;
        if vol1>1e-5
            Phase.(['fib',num2str(i)])=PhaseData(['fib',num2str(i)],Mat.fib,{'spheroid',ar,[leb.azi(i),leb.zeni(i)]},v.fib*vol1);
            Phase.(['fib',num2str(i)]).IF=comp.IF.par;
            str(j)=eval(['Phase.fib',num2str(i)]);
            j=j+1;
        else
            str(j)=[];
        end
    end
else
    error('Fiber orientation not properly defined')
end

% CONSIDER ASPECT RATIO DISTRIBUTION
nfam1=20; % number of families related to aspect ratio distribution
if strcmp(comp.fib.ardis,'single')
    % nothing to do
else
    nfam=numel(str); %number of families related to orientation distribution
    yl=linspace(0+1/nfam1/2,1-1/nfam1/2,nfam1); %equidistant probabilities
    %
    if strcmp(comp.fib.ardis,'uniform')
        a=comp.fib.ar(1);b=comp.fib.ar(2);
        xl=unifinv(yl,a,b);
    elseif strcmp(comp.fib.ardis,'weibull')
        lambda=comp.fib.ar(1);k=comp.fib.ar(2);
        xl=wblinv(yl,lambda,k);
    elseif strcmp(comp.fib.ardis,'normal')
        mu=comp.fib.ar(1);si=comp.fib.ar(2);
        xl=norminv(yl,mu,si);
    elseif strcmp(comp.fib.ardis,'lognormal')
        mu=comp.fib.ar(1);si=comp.fib.ar(2);
        xl=logninv(yl,mu,si);
    end
    fibcount=0;
    %str=repmat(str,1,nfam1);
    for i=1:nfam
        for j=1:nfam1
            if isfield(Phase,['fib',num2str(i)])
            fibcount=fibcount+1;                
            tmp=copy(Phase.(['fib',num2str(i)]));
            tmp.name=[tmp.name,'ar',num2str(j)];
            tmp.shape.sr=1/xl(j);
            tmp.vol=tmp.vol/nfam1;
            str(fibcount)=tmp;
            end
        end
    end
end

RVE.comp=RVEData1('comp',[Phase.mat,str,Phase.air],'MT',{'matrix'},tol);
[RVE.comp,Phase.comp]=homE(RVE.comp);
end

