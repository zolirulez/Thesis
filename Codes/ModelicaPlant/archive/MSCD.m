function [g,fault,t,threshold] = MSCD(x,n,threshold)
% Moving square change detector
r = NaN(length(x)-n,1);
m = movmean(x,n);
g = NaN(length(x)-n,1);
r(1) = var(x(1:n));
m(1) = mean(x(1:n));
g(1) = 0;

for it = 2:length(x)-n
%     m(it) = m(it-1) + 1/n*(x(it+n)-x(it));
    %r(it) = r(it-1) + 1/n*((x(it+n)*m(it+n))^2-(x(it)*m(it))^2); %- (m(it)^2-m(it-1)^2));
    r(it) = r(it-1) + 1/n*(x(it+n)^2-x(it)^2);
    g(it) = r(it); % max(0,g(it-1)+(sign(m(it))*m(it)-sqrt(v(it))));
end

if nargin < 3
    threshold = 1.5*max(r(1:round(length(r)/5)));
end

fault = g > threshold;
t = n+1:length(x);