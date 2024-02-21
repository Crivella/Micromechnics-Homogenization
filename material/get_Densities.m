function [rho]=get_Densities(M,Study,numcomp)

rho.C3S=3.15+zeros(1,numcomp);
rho.C2S=3.27+zeros(1,numcomp);
rho.C3A=3.03+zeros(1,numcomp);
rho.C4AF=3.708+zeros(1,numcomp);
rho.H=1.00+zeros(1,numcomp);
rho.CH=2.25+zeros(1,numcomp);
rho.ETT=1.78+zeros(1,numcomp);
rho.MSA=2.015+zeros(1,numcomp); % see NANOCEM
rho.TCA=2.042+zeros(1,numcomp); % see NANOCEM
rho.TALC=1.945+zeros(1,numcomp);
if Study==1;
    rho.SLAG=[1 2.89 1 2.89 1];
elseif Study==2;
    rho.SLAG=[1 2.93 2.91 1 2.93 2.91 1 2.93 2.91];
elseif Study==3;
    rho.SLAG=[1 1 1 2.95 2.95 2.95 2.95 2.95 2.95];
elseif Study==4;
    rho.SLAG=[1 2.93 1 1];    
end
rho.LMF=2.714+zeros(1,numcomp);
rho.LF=2.71+zeros(1,numcomp);
rho.HEMHYD=2.76+zeros(1,numcomp);
rho.GYPS=2.31+zeros(1,numcomp);
rho.ANH=2.97+zeros(1,numcomp);
rho.CEMfill=rho.HEMHYD;
rho.SLAGfill=rho.SLAG;
rho.freeM=3.58+zeros(1,numcomp);
rho.freeC=3.35+zeros(1,numcomp);
rho.Quartz=2.65+zeros(1,numcomp);
rho.FH=3+zeros(1,numcomp);
rho.sCSH=2.604+zeros(1,numcomp);
rho.sCASH=2.73+zeros(1,numcomp); % both values come from SANS/SAXS studies from Thomas,Allen,Jennings
rho.CSH=M.CSH./(M.sCSH./rho.sCSH+(M.CSH-M.sCSH)./rho.H);
rho.CASH=M.CASH./(M.sCASH./rho.sCASH+(M.CASH-M.sCASH)./rho.H);
end

