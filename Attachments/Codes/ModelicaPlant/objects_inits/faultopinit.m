% Script for copying variables to create fault considering estimator
faultOperation = 0;
Tfault = 0;
hFault = 0;
kff = copy(kf);
rlsf = copy(rls);
Xsf = Xs;
Yf = Y;
Uf = U;
uyf = uy;
hBPf = hBP;
Wf = W;
d1f = d1;
TBPf = TBP;
detectiontime = 0;
switchofftime = 0;

