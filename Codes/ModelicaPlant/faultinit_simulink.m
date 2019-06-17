% Fault Detector
fdGLR = FaultDetector;
meanFD.m0 = 0;
if exist('fielddata')
    variance = 6e3; %1.8614e+06;
else
    variance = 0.25;
end
param.WindowLength = 50;
param.InitialGuess = [-400,100];
param.FalseAlarmProbability = 1e-20;
method = 'GLR';
fdGLR.initialize(meanFD,variance,method,param);
clear param
fdCUSUM = FaultDetector;
meanFD.m1 = -1;
param.FalseAlarmTime = 1e20;
param.InitialGuess = [-400,100];
% param.Threshold = 50;
method = 'CUSUM';
fdCUSUM.initialize(meanFD,variance,method,param);
clear param
fdEM = FaultDetector;
method = 'EM';
param.MeanDensityRatio = 0.2;
param.SamplingTime = 1;
if ~exist('fielddata')
    param.ResponsibilityTimeConstant = 100;
else
    param.ResponsibilityTimeConstant = 10;
end
variance = diag([2e7; 2]); % 2e7
meanFD.m0 = zeros(size(rlsInitial.t,2),1);
fdEM.initialize(meanFD,variance,method,param);

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
switchofftime = 0;