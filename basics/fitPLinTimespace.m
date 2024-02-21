function [Ev,beta] = fitPLinTimespace(t,J,chooseplot)
%provides Ev and beta parameters for fitted Power Law in LC space
%   input: t...time
%          J...creep function in time space
%   output:Ev...creep modulus (for tref=86400)
%          beta...power law exponent
if nargin==2
    chooseplot=0;
end

% Set up fittype and options. for power law fitting
ft = fittype( 'power2' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.StartPoint = [50 0.2];

% Fit Power Law in time space
tref=86400;
[xData, yData] = prepareCurveData(t/tref,J);
[fitresult, ~] = fit( xData, yData, ft, opts );

% calculate Ec,beta
beta=fitresult.b;
Ev=1/fitresult.a;

if chooseplot==1
    %Plot
    figure(1)
    plot(t/tref,J,'kx')
    hold on
    plot(fitresult)
    text(0.8,0.1,['Ev=',num2str(Ev),'  beta=',num2str(beta)],'units','normalized');
    hold off
end
end

