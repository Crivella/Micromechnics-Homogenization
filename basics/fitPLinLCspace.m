function [Ev,beta] = fitPLinLCspace(p,J,chooseplot)
%provides Ev and beta parameters for fitted Power Law in LC space
%   input: p...frequency
%          J...creep function in LC space
%   output:Ev...creep modulus (for tref=86400)
%          beta...power law exponent
if nargin==2
    chooseplot=0;
end

% Set up fittype and options. for power law fitting
ft = fittype( 'power1' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.StartPoint = [0.00044 -0.1];

% Fit Power Law in LC space
[xData, yData] = prepareCurveData(p,J);
[fitresult, ~] = fit( xData, yData, ft, opts );

% calculate Ec,beta
beta=-fitresult.b;
tref=86400;
Ev=(1/tref)^beta*gamma(beta+1)/fitresult.a;

if chooseplot==1
    pplot=logspace(log10(min(p)),log10(max(p)),100);
    %Plot
    figure
    plot(p,J,'kx')
    hold on; grid on
    plot(pplot,fitresult(pplot),'r-')
    text(0.1,0.1,['Ev=',num2str(Ev),'  beta=',num2str(beta)],'units','normalized');
    set(gca,'xscale','log','yscale','log')
    hold off
end
end

