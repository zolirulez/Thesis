% Script for initialization of fault detectors

% ---------- GLR ---------------
fdGLR = FaultDetector;
meanFD.m0 = 0;
param.WindowLength = 100;
meanFD.m1 = 1e4;
% Variances are better bigger: residuals are not entirely white
if exist('fielddata')
    variance = 1e6; 
    % 20 percent should be faulty in the whole window to fire
    param.Threshold = 1/(2*variance*param.WindowLength)*(meanFD.m1*param.WindowLength*0.2)^2;
else
    variance = 2e6;
    % 20 percent should be faulty in the whole window to fire
    param.Threshold = 1/(2*variance*param.WindowLength)*(meanFD.m1*param.WindowLength*0.2)^2;
end
param.InitialGuess = [meanFD.m1,param.Threshold];
param.FalseAlarmProbability = 1e-4;
method = 'GLR';
fdGLR.initialize(meanFD,variance,method,param);
% ---------- CUSUM ---------------
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
param.SamplingTime = Ts;
param.WindowLength = 10;
% Here variance should be interpreted as a maximum bound
if exist('fielddata')
    param.MeanDensityRatio = 0.1;
    param.ResponsibilityTimeConstant = 100;
    variance = diag([1e6; 1e7])*2^2; 
else
    param.MeanDensityRatio = 0.1;
    param.ResponsibilityTimeConstant = 10;
    if enableBatch
        variance = diag([1e7; 1e7])*2^2; 
    else
        variance = diag([2e6; 1e7])*2^2; 
    end
end
meanFD.m0 = zeros(size(rlsInitial.t,2),1);
fdEM.initialize(meanFD,variance,method,param);

