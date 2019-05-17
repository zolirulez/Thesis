function [Pmissed,Pfalse] = GLR_design(M, mu_1,mu_0,sigma,threshold,potplease)
% script to illustrate GLR design for mean value 
% change detection of Gaussian signal
% DFTC 3rd ed. 2015, pp 291 - 293
%
% Let M be window size
if(nargin <1)
    M = 50;
end
if(nargin < 2)
    mu_1 = 0.5;
end
if(nargin < 3)
    mu_0 = 0.0;
end
if(nargin <4)
    sigma = 1.0;
end
% lambda is the non-centrality parameter in the non central 
% Chi square distribution, see (7.39).
% The GLRT has the test statistic given by (7.35)
% where z(i) is Normal with mean mu_1 and variance sigma^2.
%
% g(k) = 1/(2*sigma^2*M)* (sum(i=k-M+1,k)(z(i)-mu_0))^2
% the noncentrality parameter for the distribution of 2*g_M is (7.39):

% range for g(k)
g = 0:0.1:15;

%Under H0:
dof = 1;
pdfH0 = chi2pdf(2*g,dof);
cdfH0 = chi2cdf(2*g,dof);

%Under H1:
dof = 1;
lambda = M *(mu_1-mu_0)^2/sigma^2;
pdfH1 = ncx2pdf(2*g,dof,lambda);
cdfH1 = ncx2cdf(2*g,dof,lambda);

% choose threshold h (threshold for 2*g is 2*h !!)

h = threshold;
% compute PF and PD:
PF = 1 - chi2cdf(2*h, dof);
PD = 1 - ncx2cdf(2*h, dof,lambda);

nullh  = [0:0.1:2*h];
PFline = (1-PF)*ones(size(nullh));
PDline = (1-PD)*ones(size(nullh));

hrangep = [0:0.1:0.3];
hlinep  = 2*h*ones(size(hrangep));
hrangec = [0:0.1:1];
hlinec  = 2*h*ones(size(hrangec));

if potplease
    % plot the result
    figure,
    plot(2*g,cdfH0,'b-',2*g,cdfH1,'r',nullh,PFline,'b--',...
        nullh,PDline,'r--',hlinec,hrangec,'k--');
    xlabel('amplitude of 2*g')
    title(['cdf under H0 and H1 for GLR test, window M =',num2str(M)]);
    legend(['H0';'H1';'PF';'PD';'2h'])
    grid
    mytext = ['\mu_1 = ',num2str(mu_1,2),', \mu_0 = 0',', \sigma = ',num2str(sigma,2)];
    text(10,0.2,mytext)
    
    % pdf olots
    figure,
    plot(2*g,pdfH0,'b-',2*g,pdfH1,'r',hlinep,hrangep,'k--');
    xlabel('amplitude of 2*g')
    title(['pdf under H0 and H1 for GLR test, window M =',num2str(M)]);
    legend(['H0';'H1';'2h'])
    grid
    mytext = ['\mu_1 = ',num2str(mu_1,2),', \mu_0 = 0',', \sigma = ',num2str(sigma,2)];
    text(10,0.2,mytext)
end

% Returns
Pmissed = 1 - PD;
Pfalse = PF;
end
