function [ res ] = fun_HSbounds(M)
%HSbounds for n-phase isotropic material
% input is a matrix with n rows and 3 columns
%   column1: k1,k2,...,kn -- bulk moduli (kn>k2>k1)
%   column2: mu1,mu2,...,mun --  bulk moduli (mun>mu2>mu1)
%   column3: f1,f2,...fn -- volume fractions
% output matrix with [klow,kup;mulow,muup]


k=M(:,1);
mu=M(:,2);
f=M(:,3);
if abs(sum(f)-1)>eps
    error('volume fractions do not add up to 1')
end

al=-3./(3*k+4*mu);
A1aux=f./(1./(k-k(1))-al(1));
A1=sum(A1aux(2:end));
Anaux=f./(1./(k-k(end))-al(end));
An=sum(Anaux(1:end-1));

klow=k(1)+A1/(1+al(1)*A1);
kup=k(end)+An/(1+al(end)*An);


be=-(3*(k+2*mu))./(5*mu.*(3*k+4*mu));
B1aux=f./(1./(2*(mu-mu(1)))-be(1));
B1=sum(B1aux(2:end));
Bnaux=f./(1./(2*(mu-mu(end)))-be(end));
Bn=sum(Bnaux(1:end-1));

mulow=mu(1)+0.5*B1/(1+be(1)*B1);
muup=mu(end)+0.5*Bn/(1+be(end)*Bn);

res=[klow,kup;mulow,muup];
end

