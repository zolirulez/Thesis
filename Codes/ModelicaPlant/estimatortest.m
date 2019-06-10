clearvars
format long
% Initializing FMIKit and adding paths
addpath('C:\Users\u375749\Documents\Thesis\Codes\Linearization')
addpath('C:\Users\u375749\Documents\Thesis\Codes\MatlabPlant')
% 'Running modelcreation_simplified.m
modelcreation_simplified
% Running substitution_simplified.m
[A,B,C,D] = constantsave(A,B,C,D,Qcont);
substitution_simplified
% Running kfinit_simulink.m
kfinit_simulink

% Experiment for observation
gFunction = ginit(y);
UFunction = u;
YFunction = y;
% Experiment for actuation
XFunction = x;
load uy_sim_faulty

U = uy_sim.signals.values(:,1:nu);
Y = uy_sim.signals.values(:,nu+1:nu+ny);
ny = size(Y,2);
uy = [zeros(nu+ny,1); reshape(system.A,nx*nx,1); reshape(system.B,nx*nu,1); reshape(system.C,ny*nx,1); reshape(system.D,ny*nu,1); reshape(noise.Q,nx*nx,1)];
ABCDQ = [reshape([system.A system.B; system.C system.D],(nx+ny)*(nx+nu),1); reshape(noise.Q,nx*nx,1)];
Xs = initial.xs;
finish = 10000;
start = 2001;
delay = 200;
% Delay of fan and compressors
% U(start-1:end,4) = U(start-1-delay:end-delay,4);
% U(start-1:end,10) = U(start-1-delay:end-delay,10);
% U(start-1:end,11) = U(start-1-delay:end-delay,11);
U(:,12) = 0.2*ones(length(U),1);
% Parameter estimation
W = [sigmaValues(1); sigmaValues(2)];
sigma = 74300/3;
nw = length(W);
P = diag(W)/1e3;
W = [5000; 6000];
% Constraints (note: feedback)
TBP = CoolProp.PropsSI('T','P',Y(start,1),'H',Y(start,2),'CO2');
TBPf = TBP;
dR = CoolProp.PropsSI('D','P',Y(start,3),'H',Y(start,4),'CO2');
d1 = CoolProp.PropsSI('D','P',Y(start,1),'H',Xs(2),'CO2');
% RLS
rlsInitial.t = [5000; 6000; 1e3; 1e4; 1e3]*[1 1e-3];  % 10
rlsInitial.P = diag(max(rlsInitial.t')')/10; 
rls = RecursiveLeastSquares;
lambda = 1 - 1e-5;
weight = eye(2);
rls.initialize(lambda,weight,rlsInitial);
rlsf = copy(rls);
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
param.MeanDensityRatio = 0.005;
param.SamplingTime = 1;
param.ResponsibilityTimeConstant = 100;
variance = diag([4.25e6; 0.25]);
mean.m0 = zeros(2,1);
fdEM.initialize(mean,variance,method,param);
% Boolean for fault operation
faultOperation = 0;
% Low pass filter
hBP = Y(start,2);
TA0 = U(start,10);
h1 = Y(start,5);
dBP = U(start,6);
Tfault = 0;
% Recording
record = NaN(nx+nx+nx+ny+nw,finish-start+1);
resrecord = NaN(2,finish-start+1);
paramrecord = NaN(size(rls.t,1)*size(rls.t,2),finish-start+1);
outrecord = NaN(2,finish-start+1);
outcorrecord = [];
grecord = NaN(3,finish-start+1);
faultrecord = NaN(3,finish-start+1);
% Fault operation Kalman Filter
kff = copy(kf);
Xsf = Xs;
Yf = Y;
uyf = uy;
ABCDQf = ABCDQ;
hBPf = Yf(start,2);
h1f = Yf(start,5);
d1f = d1;
detectiontime = 0;
recordf = record;
for it = start:finish
    if ~rem(it,100)
        disp(['Iteration ' num2str(it)])
    end
    % ------------ State estimation -----------
    xf = kf.markovPredictor(uy(1:nu),uy(nu+1:nu+ny),reshape(uy(nu+ny+1:nu+ny+nx*nx),nx,nx),reshape(uy(nu+ny+nx*nx+1:nu+ny+nx*nx+nx*nu),nx,nu),reshape(uy(nu+ny+nx*nx+nx*nu+1:nu+ny+nx*nx+nx*nu+ny*nx),ny,nx), reshape(uy(nu+ny+nx*nx+nx*nu+ny*nx+1:nu+ny+nx*nx+nx*nu+ny*nx+ny*nu),ny,nu), reshape(uy(nu+ny+nx*nx+nx*nu+ny*nx+ny*nu+1:nu+ny+nx*nx+nx*nu+ny*nx+ny*nu+nx*nx),nx,nx));
    xff = kff.markovPredictor(uyf(1:nu),uyf(nu+1:nu+ny),reshape(uyf(nu+ny+1:nu+ny+nx*nx),nx,nx),reshape(uyf(nu+ny+nx*nx+1:nu+ny+nx*nx+nx*nu),nx,nu),reshape(uyf(nu+ny+nx*nx+nx*nu+1:nu+ny+nx*nx+nx*nu+ny*nx),ny,nx), reshape(uyf(nu+ny+nx*nx+nx*nu+ny*nx+1:nu+ny+nx*nx+nx*nu+ny*nx+ny*nu),ny,nu), reshape(uyf(nu+ny+nx*nx+nx*nu+ny*nx+ny*nu+1:nu+ny+nx*nx+nx*nu+ny*nx+ny*nu+nx*nx),nx,nx));
    Xs = kf.x1 + Xs;
    Xsf = kff.x1 + Xsf;
    Xs(3) = d1; % TODO
    Xs(7) = dR; % TODO
    Xsf(3) = d1f; % TODO
    Xsf(7) = dR; % TODO
    % ------------ Parameter estimation -----------
    % Delayed THR!
    THR = CoolProp.PropsSI('T','H',U(it-delay,11),'P',Y(it-delay,1),'CO2');
    CRIT = U(it-1,4);
    CRA = U(it-1,1);
    dBP = U(it+1,6);
    DmV = U(it,3)*KvValues*sqrt(dBP*(Y(it,1) - Y(it,3)));
    % Delayed hHR!
    DQ = DmV*(U(it-delay,11) - Y(it,2));
    DQf = DmV*(U(it-delay,11) - Yf(it,2));
    TA0 = U(it+1,10);
    DT = max(1,TBP - TA0);
    DTf = max(1,TBPf - TA0);
    phi = [1; CRA; TA0; CRIT; THR];
    out = [DQ DT]; % Xs(8)
    outf = [DQf DTf];
    rls.regression(phi,out);
    rlsf.regression(phi,outf);
    W = [DQ-phi(2)*rls.t(2)'; rls.t(2)]/DT;
    Wf = [DQf-phi(2)*rlsf.t(2)'; rlsf.t(2)]/DTf;
    % ------------ Fault Operation -----------
    if it == start
        rls.e = zeros(1,length(rls.e));
        [~,fault1] = fdCUSUM.CUSUM(rls.e(2)); % TODO: wrong place
        [~,fault2] = fdGLR.GLR(rls.e(2)); % TODO: wrong place
    end
    [~,fault3] = fdEM.EM(rls.e');
    % Measurement correction: based on CUSUM
    if it > start + 20
        if all(faultrecord(3,it-19-start:it-start)) && ~faultOperation
            Wsave = reshape(paramrecord(:,it-300-start),length(phi),length(out));
            if ~detectiontime
                detectiontime = it;
            end
            faultOperation = 1;
        end
        if ~any(faultrecord(3,it-19-start:it-start))
            faultOperation = 0;
        end
    end
    if faultOperation
        outcor = phi'*Wsave;
        outcorrecord = [outcorrecord; outcor];
        Tfault = out(2) - outcor(2);
    else
        Tfault = 0;
    end
    % ------------ Parameter substitutions for new iteration -----------
    UXYW = [U(it,:)'; Xs; Y(it,1:ny)'; W];
    UXYWf = [U(it,:)'; Xsf; Yf(it,1:ny)'; Wf];
    ABCDQ = LTVsystemDescription(UXYW(1:nu), UXYW(nu+1:nu+nx), UXYW(nu+nx+1:nu+nx+ny), UXYW(nu+nx+ny+1:nu+nx+ny+nw));
    ABCDQf = LTVsystemDescription(UXYWf(1:nu), UXYWf(nu+1:nu+nx), UXYWf(nu+nx+1:nu+nx+ny), UXYWf(nu+nx+ny+1:nu+nx+ny+nw));
    A = ABCDQ(1:nx*nx);
    B = ABCDQ(nx*nx+1:nx*nx+nx*nu);
    C = ABCDQ(nx*nx+nx*nu+1:nx*nx+nx*nu+ny*nx);
    D = ABCDQ(nx*nx+nx*nu+ny*nx+1:nx*nx+nx*nu+ny*nx+ny*nu);
    Q = ABCDQ(nx*nx+nx*nu+ny*nx+ny*nu+1:nx*nx+nx*nu+ny*nx+ny*nu+nx*nx);
    Af = ABCDQf(1:nx*nx);
    Bf = ABCDQf(nx*nx+1:nx*nx+nx*nu);
    Cf = ABCDQf(nx*nx+nx*nu+1:nx*nx+nx*nu+ny*nx);
    Df = ABCDQf(nx*nx+nx*nu+ny*nx+1:nx*nx+nx*nu+ny*nx+ny*nu);
    Qf = ABCDQf(nx*nx+nx*nu+ny*nx+ny*nu+1:nx*nx+nx*nu+ny*nx+ny*nu+nx*nx);
    % ------------ Constraints -----------
    TBP = CoolProp.PropsSI('T','P',Y(it+1,1),'H',Y(it+1,2),'CO2');
    dR = CoolProp.PropsSI('D','P',Y(it+1,3),'H',Y(it+1,4),'CO2');
    d1 = CoolProp.PropsSI('D','P',Y(it+1,1),'H',Y(it+1,5),'CO2');
    % In case of a fault
    TBPf = TBP - Tfault;
    Yf(it+1,2) = CoolProp.PropsSI('H','P',Yf(it+1,1),'T',TBPf,'CO2');
    d1f = CoolProp.PropsSI('D','P',Yf(it+1,1),'H',Yf(it+1,5),'CO2');
    % Low pass filters for noisy measurements
%     TauNoise = 1;
%     hBP = hBP*(1-Ts/TauNoise) + Ts/TauNoise*Y(it+1,2);
%     hBPf = hBPf*(1-Ts/TauNoise) + Ts/TauNoise*Yf(it+1,2);
%     h1 = h1*(1-Ts/TauNoise) + Ts/TauNoise*Y(it+1,5);
%     h1f = h1f*(1-Ts/TauNoise) + Ts/TauNoise*Yf(it+1,5);
%     Y(it+1,2) = hBP;
%     Yf(it+1,2) = hBPf;
%     TA0 = U(it+1,10);
%     Y(it+1,5) = T1 - (TBP - TA0);
%     Yf(it+1,5) = T1f - (TBPf - TA0);
    dBP = U(it+1,6);
    dBPf = CoolProp.PropsSI('D','P',Yf(it+1,1),'H',Yf(it+1,2),'CO2');
    U(it+1,6) = dBPf;
    % ------------ New setpoint for linearization -----------
    u = U(it+1,:)' - U(it,:)';
    y = Y(it+1,1:ny)' - g(gFunction,XFunction,Xs,UFunction,U(it,:));
    yf = Yf(it+1,1:ny)' - g(gFunction,XFunction,Xsf,UFunction,U(it,:));
    kf.x1 = zeros(nx,1);
    kff.x1 = zeros(nx,1);
    uy = [u; y; A; B; C; D; Q];
    uyf = [u; yf; Af; Bf; Cf; Df; Qf];
    % ------------ Recording -----------
    statecorrection = kf.Kx*kf.e;
    statecorrectionf = kff.Kx*kff.e;
    eigP1 = zeros(nx,1);%eig(kf.P1); TODO
    eigP1f = zeros(nx,1);%eig(kff.P1); TODO
    record(:,it-start+1) = [Xs; eigP1; statecorrection; kf.e; W];
    recordf(:,it-start+1) = [Xsf; eigP1f; statecorrectionf; kff.e; Wf];
    resrecord(:,it-start+1) = rls.e;
    paramrecord(:,it-start+1) = reshape(rls.t,size(rls.t,1)*size(rls.t,2),1);
    outrecord(:,it-start+1) = out';
    grecord(:,it-start+1) = [fdCUSUM.g; fdGLR.g; fdEM.g];
    faultrecord(:,it-start+1) = [fault1; fault2; fault3];
end

plotting