function sketch4em(resid,LR,var0,tau,Ts)
% usage: sketch4em(testresid,0.999,var(resid_normal),100,1)
% tau: responsibility dynamics
% Ts: sampling time
% var0: variance of normal residual
% LR: likelihood ratio of fault at fault distribution mean value
m0 = 0;
m1 = m0 + sqrt(-2*var0*log(1/LR - 1));
pi0 = 0.95;
pi1 = 1-pi0;
pi1v = zeros(length(resid),1);
vv = zeros(length(resid),1);
var1 = var0;
m = 0;
for it = 2:length(resid)
    z = resid(it);
    gain0 = (pi0+0.01)/sqrt(2*pi*var0);
    gain1 = (pi1+0.01)/sqrt(2*pi*var1);
    r1 = gain1*exp(-0.5*(z-m1)^2/var1);
    r0 = gain0*exp(-0.5*(z-m0)^2/var0);
%     pi1 = r1/(r0 + r1);
pi1 = (1-Ts/tau)*pi1 + Ts/tau*r1/(r0 + r1);
    pi0 = 1-pi1;
    mn = pi1*m1 + pi0*m0;
%     m = tau/Ts*(mn - (1-Ts/tau)*m);
%     pi0 =  (m - m1)/(m0 - m1);
%     pi1 = 1-pi0;
    var1 = (z-mn)'*(z-mn);
    pi1v(it) = pi1;
    vv(it) = var1;
end
figure, plot(pi1v)
figure, plot(pi1v>0.5)