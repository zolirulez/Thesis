% Fault Detector
fdGLR = FaultDetector;
mean.m0 = 0;
variance = 0.25; %1.8614e+06;
param.WindowLength = 500;
param.InitialGuess = [-1,1];
param.FalseAlarmProbability = 1e-20;
method = 'GLR';
fdGLR.initialize(mean,variance,method,param);
clear param
fdCUSUM = FaultDetector;
mean.m1 = -1;
param.FalseAlarmTime = 1e20;
param.InitialGuess = [-1,1];
method = 'CUSUM';
fdCUSUM.initialize(mean,variance,method,param);
clear param
fdEM = FaultDetector;
method = 'EM';
param.MeanDensityRatio = 0.2;
param.SamplingTime = 1;
param.ResponsibilityTimeConstant = 100;
variance = diag([2e7; 2]); % 2e7
mean.m0 = zeros(size(rlsInitial.t,2),1);
fdEM.initialize(mean,variance,method,param);

% Fault operation
faultOperation = 0;
Tfault = 0;
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