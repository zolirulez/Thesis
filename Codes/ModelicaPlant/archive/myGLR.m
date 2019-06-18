function [g,fault,threshold] = myGLR(mu0,var,resid,M,threshold)
% The function is given by
%           g(k) = 1/(2*sigma^2*M)* (sum(i=k-M+1,k)(z(i)-mu_0))^2
%
% g is the decision function
% Var is the variance of the residual
% The residual is a timeseries vector
% IMPORTANT: if threshold is not provided, then the first half of
%       the timeseries is to be in a faultless phase

if isstruct(resid)
    resid = resid.Data';
end
r = [mu0*ones(M,1); resid'];
g = zeros(length(resid')+M,1);
for k = M+1:length(r)
    g(k) = 1/(2*var*M)*(sum(r(k-M+1:k)-mu0))^2;
end
g = g(M+1:end);
if nargin<5
    threshold = 1.1*max(abs(g(1:round(length(g)/5))));
end
fault = (g>threshold);