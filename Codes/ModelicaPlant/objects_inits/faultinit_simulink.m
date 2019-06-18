% Script for initialization of fault detectors

% ---------- GLR ---------------
fdGLR = FaultDetector;
meanFD.m0 = 0;
param.WindowLength = 100;
% Variances are better bigger, like 5 times: residuals are not entirely
% white
if exist('fielddata')
    variance = 5e6; % 6e3; %1.8614e+06;
    param.InitialGuess = [-sqrt(variance)*5,param.WindowLength/variance];
    param.Threshold = param.WindowLength/variance*1000;
else
    variance = 5e6;
    param.InitialGuess = [sqrt(variance)*5,param.WindowLength/variance];
    param.Threshold = param.WindowLength/variance*5;
end
param.FalseAlarmProbability = 1e-4;
method = 'GLR';
fdGLR.initialize(meanFD,variance,method,param);
% ---------- CUSUM ---------------
% CUSUM does not work currently
clear param
fdCUSUM = FaultDetector;
if exist('fielddata')
    meanFD.m1 = -sqrt(variance)*5;
else
    meanFD.m1 = sqrt(variance)*5;
end
% param.FalseAlarmTime = 1e20;
% param.InitialGuess = [-400,100];
if exist('fielddata')
    param.Threshold = 1/variance;
else
    param.Threshold = 1/variance;
end
method = 'CUSUM';
fdCUSUM.initialize(meanFD,variance,method,param);
% ---------- EM ---------------
clear param
fdEM = FaultDetector;
method = 'EM';
param.MeanDensityRatio = 0.2;
param.SamplingTime = Ts;
% Here variance should rather be interpreted as a maximum bound
if exist('fielddata')
    param.ResponsibilityTimeConstant = 10;
    variance = diag([2e7; 2]); 
else
    param.ResponsibilityTimeConstant = 10;
    variance = diag([2e7; 2]); % 2e7
end
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