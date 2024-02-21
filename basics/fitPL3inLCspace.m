function [fitr] = fitPL3inLCspace(p,J,chooseplot)
%provides Ev and beta parameters for fitted Power Law in LC space
%   input: p...frequency
%          J...creep function in LC space
%   output:Ev...creep modulus (for tref=86400)
%          beta...power law exponent
if nargin==2
    chooseplot=0;
end

tref=86400;
% Set up fittype and options. for power law fitting
ft = fittype( 'power2' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.StartPoint = [0.00044 -0.1 0];

% Fit Power Law in LC space
[xData, yData] = prepareCurveData(p,J);
[fitresult, ~] = fit( xData, yData, ft, opts );

% calculate parameter for time function f(t)=a*t^b+c
fitr.b=-fitresult.b;
fitr.a=(1/tref)^fitr.b*fitresult.a/gamma(fitr.b+1);
fitr.c=fitresult.c;

if chooseplot==1
    pplot=logspace(log10(min(p)),log10(max(p)),100);
    %Plot
    figure
    plot(p,J,'kx')
    hold on; grid on;
    plot(pplot,fitresult(pplot),'r-')
    text(0.1,0.1,['a=',num2str(fitr.a),'  b=',num2str(fitr.b),'  c=',num2str(fitr.c)],'units','normalized');
    set(gca,'xscale','log','yscale','log')
    hold off
end
end
