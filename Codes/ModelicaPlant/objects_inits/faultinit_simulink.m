% Script for initialization of fault detectors

% ---------- GLR ---------------
fdGLR = FaultDetector;
meanFD.m0 = 0;
param.WindowLength = 100;
if exist('fielddata')
    meanFD.m1 = -250; %-sqrt(variance)*5;
else
    meanFD.m1 = 2e4; % sqrt(variance)*5;
end
% Variances are better bigger: residuals are not entirely white
if exist('fielddata')
    variance = 4.5e5; % 6e3; %1.8614e+06;
    % 20 percent should be faulty in the whole window to fire
    param.Threshold = 1/(2*variance*param.WindowLength)*(meanFD.m1*param.WindowLength*0.2)^2;
else
    variance = 5e6;
    param.Threshold = 1/(2*variance*param.WindowLength)*(meanFD.m1*param.WindowLength*0.2)^2;
end
param.InitialGuess = [meanFD.m1,param.Threshold];
param.FalseAlarmProbability = 1e-4;
method = 'GLR';
fdGLR.initialize(meanFD,variance,method,param);
% ---------- CUSUM ---------------
% CUSUM does not work currently
clear param
fdCUSUM = FaultDetector;
% 20 pieces of mean value should be measured in order to fire
param.FalseAlarmTime = 1e4;
param.Threshold = meanFD.m1^2/2/variance*20;
param.InitialGuess = [meanFD.m1,param.Threshold];
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