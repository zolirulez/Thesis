FMIKit.initialize
addpath('C:\Users\u375749\Documents\Thesis\Codes\Linearization')
addpath('C:\Users\u375749\Documents\Thesis\Codes\MatlabPlant')
Ts = 1;
negation = -1;
minOutput = 0;
maxOutput = 1;
PI_HP = PIController;
Ti = 50;
K = 1e-7;
initialIntegration = 0.25;
PI_HP.initialize(K,Ti/Ts,initialIntegration,minOutput,maxOutput,negation,Ts)
%  K_Controllers.K_PI k_PI_HighPressure(
%     Ti=50,
%     yMin=0,
%     k=1e-5,
%     yInitial=0.25)
%     annotation (Placement(transformation(extent={{6,-6},{-6,6}},
%         rotation=90,
%         origin={-112,156})));
PI_GC = PIController;
Ti = 200;
K = 1e-2;
initialIntegration = 0.625;
PI_GC.initialize(K,Ti/Ts,initialIntegration,minOutput,maxOutput,negation,Ts)
%   K_Controllers.K_PI k_PI_GCOutlet(
%     Ti=200,
%     yInitial=0.625,
%     k=1e-2)
PI_MT = PIController;
Ti = 100;
K = 1e-6;
initialIntegration = 0.2;
PI_MT.initialize(K,Ti/Ts,initialIntegration,minOutput,maxOutput,negation,Ts)
%   K_Controllers.K_PI k_PI_MTCompressor(
%     Ti=100,
%     k=1e-6,
%     yInitial=0.2)
PI_LT = PIController;
Ti = 50;
K = 1e-7;
initialIntegration = 0.3;
PI_LT.initialize(K,Ti/Ts,initialIntegration,minOutput,maxOutput,negation,Ts)
%   K_Controllers.K_PI k_PI_LTCompressor(
%     Ti=50,
%     k=1e-7,
%     yInitial=0.3)
PI_RP = PIController;
Ti = 20;
K = 1e-5;
initialIntegration = 0.3;
PI_RP.initialize(K,Ti/Ts,initialIntegration,minOutput,maxOutput,negation,Ts)
%   K_Controllers.K_PI k_PI_ReceiverPressure(Ti=20,
%     k=1e-4,
%     yInitial=0.3)
PI_BP = PIController;
Ti = 10;
K = 1e-4;
initialIntegration = 1;
PI_BP.initialize(K,Ti/Ts,initialIntegration,minOutput,maxOutput,negation,Ts)
%   K_Controllers.K_PI k_PI_ByPass(
%     Ti=10,
%     yInitial=1,
%     k=1e-4)
PI_C = PIController;
Ti = 100;
K = 1e-3;
initialIntegration = 0.1;
PI_C.initialize(K,Ti/Ts,initialIntegration,minOutput,maxOutput,negation,Ts)
%   K_Controllers.K_PI k_PI_Cooler(
%     Ti=100,
%     k=1e-3,
%     yInitial=0.1)
PI_F = PIController;
Ti = 100;
K = 1e-3;
initialIntegration = 0.04;
PI_F.initialize(K,Ti/Ts,initialIntegration,minOutput,maxOutput,negation,Ts)
%   K_Controllers.K_PI k_PI_Freezer(
%     Ti=100,
%     k=1e-3,
%     yInitial=0.04)
