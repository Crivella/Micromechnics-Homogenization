%% BIOCOMPOSITE STIFFNESS HOMOGENIZATION
% Update and cleaned up, Februar 2024

addpath(genpath('material'))
addpath(genpath('homogenization'))
addpath(genpath('basics'))
addpath(genpath('hill'))

clear; clc; close all;
warning off %otherwise constant Mori-Tanaka issue

%% Input
%note: if nothing specified, it is always a single scalar
%A)numerical
tol=8; %6 for quick computation, better 8 or even 10

%B)composite
comp.name='testcomp1'; %Name, string
comp.case=3; %either 1, 2, or 3: 1 is refers to minimum and 2 to maximum stiffness properties, 3 should be a representive average

%C)fiber macroscopic
comp.fib.name='Flax'; %mandatory, name of the fiber; to find it in the database
comp.fib.rho=NaN; %optional, fiber density [g/cm³]
%comp.fib.E=NaN; %optional, fiber modulus from experiments to compare [GPa]
%comp.fib.ft=NaN; %optional, fiber tensile strength from experiment to compare [MPa]
comp.fib.l=NaN; %optional, average fiber length [mm]
comp.fib.d=NaN; %optional, average fiber diameter [µm]
comp.fib.ardis='single'; %optional, aspect ratio distribution: "single","normal","lognormal", "weibull", "uniform", default is single
comp.fib.ar=[100];%2-element vector, optional, defining the aspect ratio distribution function, could be [ar, NaN] with ar as the mean aspect ratio
comp.fib.ori='vMs'; %mandatory, either '1D, '3Diso', '2Diso', 'vMF', "vMs"
% comp.fib.ori='align'; %mandatory, either '1D, '3Diso', '2Diso', 'vMF', "vMs"
comp.fib.oripar=[10]; %orientation distribution parameter (concentration), kappa [-]
% comp.fib.oripar=[1,1]; %orientation distribution parameter (concentration), kappa [-]

%D)fiber microscopic
comp.fib.MFA=0; %optional, microfibril angle [degrees]
comp.fib.lumpor=NaN; %optional, lumen porosity [-]
comp.fib.mfrac=NaN; %optional/mandatory, mass fraction of fiber [-], note that either mass or volume fraction has to be given
comp.fib.vfrac=1; %optional/mandatory, volume fraction of fiber [-], note that either mass or volume fraction has to be given
comp.fib.mchem.totcel=NaN; %optional, total cellulose (amorph+crystalline) content [-], w.r.t. cell wall mass, note if one is given, all have to be given
comp.fib.mchem.hemcel=NaN; %optional, hemicellulose content [-], w.r.t. cell wall mass
comp.fib.mchem.lignin=NaN; %optional, lignin content [-], w.r.t. cell wall mass
comp.fib.mchem.pectin=NaN; %optional, pectin content [-], w.r.t. cell wall mass
comp.fib.mchem.extr=NaN; %optional, extractives content [-], w.r.t. cell wall mass
comp.fib.mchem.wax=0; %optional, wax content [-], w.r.t. cell wall mass
comp.fib.mchem.ash=NaN; %optional, ash content [-], w.r.t. cell wall mass
comp.fib.CI=NaN; %optional, crystallinity of cellulose [-]

%E)matrix
comp.mat.name='PP';%mandatory, name of the matrix; to find it in the database
comp.mat.rho=NaN;%optional, matrix density [g/cm³], if fiber volume fraction not given
comp.mat.E=5;%optional, matrix modulus [GPa]
comp.mat.nu=-0,8;%optional, matrix Poisson's ratio [-], note that stiffness of matrix has to be defined, so if not in the database, it has to be defined here
%comp.mat.ft=NaN;%optional, matrix tensile strength from experiment to compare [MPa]
%comp.mat.eps=NaN;%optional, matrix ultimate strain from experiment to compare [MPa]
comp.mat.fair=0.05;%optional, air porosity [-], default is zero


%F) interface
comp.IF.mode=0; % empty or 0...perfect bond, 1...imperfect bond
comp.IF.par=[0 50]; % two-element vector with interface compliances [alpha, beta], tangential and longitudinal compliance in 1/GPa, attention: normalization w.r.t. fiber radius, actually m/GPa if not normalized in R function

% additional parameters
comp.MLar=1e-20; % aspect ratio of Middle lamella/lumen, either infinitely long 1e-20, or smaller [-]

%% define composition/microstructure parameters based on input
comp=prep_comp_noexp(comp);


%% Homogenization
[RVE]=hom_biocomp_stiffness(comp,tol);
comp.hom.C=RVE.comp.Chom;
comp.hom.Cfib=RVE.fib.Chom;
Dhom=inv(RVE.fib.Chom);
comp.hom.Efib=1/Dhom(3,3);
if strcmp('iso',fun_check_isotropy(RVE.comp.Chom))
    [comp.hom.El,comp.hom.nu]=fun_Enu_from_C(RVE.comp.Chom);
    comp.hom.Et=comp.hom.El;
else
    Dhom=inv(comp.hom.C);
    if strcmp(comp.fib.ori,'2Diso') || strcmp(comp.fib.ori,'vMs') ;
        comp.hom.El=1/Dhom(1,1);comp.hom.nu=NaN;
        comp.hom.Et=1/Dhom(3,3);
    elseif strcmp(comp.fib.ori,'vMF') || strcmp(comp.fib.ori,'align') || strcmp(comp.fib.ori,'1D') || strcmp(comp.fib.ori,'aligned');
        comp.hom.El=1/Dhom(3,3);comp.hom.nu=NaN;
        comp.hom.Et=1/Dhom(1,1);
    else error()
    end
end
disp(['Comp:',comp.name,' c',num2str(comp.case),': El=',num2str(comp.hom.El,3),', Et=',num2str(comp.hom.Et,3)])
