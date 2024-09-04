Mat.PP=MatData('PP','iso',{'Enu',1.18, 0.43},false,false); 
Mat.PLA=MatData('PLA','iso',{'Enu',100, 0.33},false,false);

Phase.PP=PhaseData('PP',Mat.PP,{'spheroid',1e-20,[0,45*pi/180]},0.6);
Phase.PLA=PhaseData('PLA',Mat.PLA,{'matrix'},0.4);

RVE=RVEData1('test',[Phase.PP,Phase.PLA],'MT',{'matrix'},8);
[RVE,Phase.hom]=homE(RVE);
