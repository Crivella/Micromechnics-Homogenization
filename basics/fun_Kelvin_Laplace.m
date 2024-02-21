function [pc,pa] = fun_Kelvin_Laplace(RH,varargin)
%fun_Kelvin_Laplace computes the capillary pore pressure pc

%default values
R=8.314; %J/mol K, ideal Gas constant
M=0.01802; %kg/mol, molar weight of water
T=293.15; %K, absolute temperature
rho=1000; %kg/m3, density of water
RHs=98; %maybe 98% or evolving?

if numel(varargin)==1
    RHs=varargin{1};
end

% pressure in capillary pores
pc=-log(RH/RHs)*rho*R*T/M/1e9; %pressure in GPa
pc(pc<0)=0;


end

