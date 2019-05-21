function [g,fault,threshold] = myCusum(mu1,mu0,var,resid,threshold)
% CUSUM algorithm, for both scalar and vector case
%
% g is the decision function
% Var is the variance of the residual
% The residual is a timeseries vector
% IMPORTANT: if threshold is not provided, then the first half of
%       the timeseries is to be in a faultless phase

if isstruct(resid)
    z = resid.Data';
else
    z = resid;
end
g = zeros(length(z),1);
g(1) = max(0,(mu1-mu0)'/var*(z(:,1)-(mu0 + mu1)/2));
for k = 2:length(z)
    g(k) = max(0,g(k-1)+(mu1-mu0)'/var*(z(:,k)-(mu0 +mu1)/2));
end
if nargin<5
    threshold = 1.1*max(abs(g(1:round(length(g)/5))));
end
fault = g > threshold;