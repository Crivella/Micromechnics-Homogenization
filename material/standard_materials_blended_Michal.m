% clinker/SCM/additions
Mat.C3S=MatData('C3S','iso',{'Enu',135, 0.3},false,false); %from: Velez/Acker, nanoindentaion
Mat.C2S=MatData('C2S','iso',{'Enu',130, 0.3},false,false); %from: Velez/Acker, nanoindentaion
Mat.C3A=MatData('C3A','iso',{'Enu',145, 0.3},false,false); %from: Velez/Acker, nanoindentaion
Mat.C4AF=MatData('C4AF','iso',{'Enu',125, 0.3},false,false); %from: Velez/Acker, nanoindentaion
Mat.SLAG=MatData('slag','iso',{'Enu',78, 0.3},false,false); %from: Nemecek_CCC_2011
Mat.FA=MatData('flyash','iso',{'Enu',105, 0.2},false,false); %from:  Smilauer_JMS_2011
Mat.LF=MatData('limestone','iso',{'Enu',76.6,0.31},false,false); %E from: Presser_JMS_2010, nu estimated
Mat.ANH=MatData('gypsum','iso',{'Enu',73.6,0.31},false,false); %from: C.-J. Haecker_CCR_2005, note that E=45.7 and nu=0.33 reported in article do not correspond to the used k and mu!
Mat.QUARTZ=MatData('quartz','iso',{'Enu',72.8,0.167},false,false); %from: Ahrens_Mineral Physics_1995

%aggregates
Mat.sand=MatData('sand','iso',{'Enu',70,0.17},false,false); %from: Vorel_JCAM_2012
Mat.agg=MatData('agg','iso',{'Enu',70,0.17},false,false); % 

% hydrates
Mat.hyd=MatData('hyd','iso',{'Enu',29.15786664, 0.24},true,{'DP',[0.06068,0.2580]}); %from: Pichler/Hellmich_CCR_2011, %DP parameter from nanoindenation+limit analysis of Sarris/Constantinides, and from calibration based on MC-paramters
Mat.CH=MatData('CH','iso',{'Enu',42.3, 0.324},false,false); %from: Haecker_CCR_2005 and older refs (Holuj_SSC_1985, Monteiro_CCR_1995)
Mat.AFt=MatData('AFt','iso',{'Enu',25.41,0.344},false,false);
Mat.AFm=MatData('AFm','iso',{'Enu',42.3, 0.324},false,false); % equivalent to CH, as done from Haecker_CCR_2005
Mat.MCA=MatData('MCA','iso',{'Enu',42.3, 0.324},false,false);
Mat.MH=MatData('MH','iso',{'Enu',42.3, 0.324},false,false);
Mat.otherhyd=MatData('otherhyd','iso',{'Enu',42.3, 0.324},false,false);

% other
Mat.pore=MatData('pore','iso',zeros(6),false,false);