% Script for testing the estimator on simulation or field data
% Zoltan Mark Pinter, Master Thesis, 2019

% clearvars
format long
addpath('C:\Users\u375749\Documents\Thesis\Codes\ModelicaPlant\bestsimulations')
addpath('C:\Users\u375749\Documents\Thesis\Codes\ModelicaPlant\objects_inits')
addpath('C:\Users\u375749\Documents\Thesis\Codes\Linearization')
modelcreation_simplified
eval(['load ' simulationstring]) % Some of modelcreation_simplified may be overwritten
if exist('fielddata')
    clear fielddata
end

% Output function
YFunction = matlabFunction(ginit(y));
UFunction = u;
XFunctionSym = x;

% Time parameters, data sets
U = uy_sim.signals.values(:,1:nu);
Y = uy_sim.signals.values(:,nu+1:nu+ny);
start = 1001;
finish = length(Y)-1;
U(U(:,5)>0.5,5) = 0.5;
hRatio = 1/3;
Y(:,5) = hRatio*[U(round(delayh/2),10)*ones(delayh,1); U(1:end-delayh,10)] + (1-hRatio)*Y(:,2);
ny = size(Y,2);

% Running initialization functions for tools
uy = [zeros(nu+ny,1); reshape(system.A,nx*nx,1); reshape(system.B,nx*nu,1); reshape(system.C,ny*nx,1); reshape(system.D,ny*nu,1); reshape(noise.Q,nx*nx,1)];

% Initial numbers
TBP = CoolProp.PropsSI('T','P',Y(start,1),'H',Y(start,2),'CO2');
dR = CoolProp.PropsSI('D','P',Y(start,3),'H',Y(start,4),'CO2');
d1 = CoolProp.PropsSI('D','P',Y(start,1),'H',Y(start,5),'CO2');
TA0 = U(start,12);
hBP = Y(start,2);
% Parameter estimation (not used now)
nw = 2;
W = [5000; 6000];

% Fault settings
faultinit_simulink
faultopinit

% Recording
record = NaN(nx+nx+ny+nw+1,finish-start+1);
resrecord = NaN(size(rlsInitial.t,2),finish+1);
paramrecord = NaN(size(rls.t,1)*size(rls.t,2),finish-start+1);
outrecord = NaN(size(rlsInitial.t,2),finish+1);
outcorrecord = NaN(size(rlsInitial.t,2)+2,finish+1);
grecord = NaN(3,finish+1);
faultrecord = NaN(3,finish+1);
recordf = record;

% Further variables
meas = Y(start,:);
measf = Yf(start,:);
constrainedEstimator = 1;
inputSample = U(start,:)';
inputSamplef = Uf(start,:)';
% Noise control (not used finally)
rng(12345)

% ------------ Fault Handling -----------
outrecord(:,:) = squeeze(yhat_sim.signals.values + residual_sim.signals.values);
outcorrecord(1:2,:) = squeeze(yhat_sim.signals.values);
outcorrecord(4,:) = h1_sim.signals.values(:,end-9);
resrecord(:,:) = squeeze(residual_sim.signals.values);
grecord(3,:) = g_sim.signals.values';
faultrecord(3,:) = g_sim.signals.values' > 0.5;
detectiontime = find(diff(faultop_sim.signals.values')==1);
try
    detectiontime = detectiontime(1);
catch
    detectiontime = find(faultrecord(3,2000:end)>0);
    detectiontime = detectiontime(1) + 2000;
end
tic
for it = start:finish
    if ~rem(it,100)
        disp(['Iteration ' num2str(it)])
    end
    % ------------ State estimation -----------
    xfsave = kf.xf;
    xffsave = kff.xf;
    kf.markovPredictor(uy(1:nu),uy(nu+1:nu+ny),reshape(uy(nu+ny+1:nu+ny+nx*nx),nx,nx),reshape(uy(nu+ny+nx*nx+1:nu+ny+nx*nx+nx*nu),nx,nu),reshape(uy(nu+ny+nx*nx+nx*nu+1:nu+ny+nx*nx+nx*nu+ny*nx),ny,nx), reshape(uy(nu+ny+nx*nx+nx*nu+ny*nx+1:nu+ny+nx*nx+nx*nu+ny*nx+ny*nu),ny,nu), reshape(uy(nu+ny+nx*nx+nx*nu+ny*nx+ny*nu+1:nu+ny+nx*nx+nx*nu+ny*nx+ny*nu+nx*nx),nx,nx));
    kff.markovPredictor(uyf(1:nu),uyf(nu+1:nu+ny),reshape(uyf(nu+ny+1:nu+ny+nx*nx),nx,nx),reshape(uyf(nu+ny+nx*nx+1:nu+ny+nx*nx+nx*nu),nx,nu),reshape(uyf(nu+ny+nx*nx+nx*nu+1:nu+ny+nx*nx+nx*nu+ny*nx),ny,nx), reshape(uyf(nu+ny+nx*nx+nx*nu+ny*nx+1:nu+ny+nx*nx+nx*nu+ny*nx+ny*nu),ny,nu), reshape(uyf(nu+ny+nx*nx+nx*nu+ny*nx+ny*nu+1:nu+ny+nx*nx+nx*nu+ny*nx+ny*nu+nx*nx),nx,nx));
    % Estimator lability handling
    if any(abs(eig(kff.A - kff.Kx*kff.C)) >= 1) 
        warning('Estimator is not stable')
        kf.x1 = xfsave;
        kff.x1 = xffsave;
    end
    % Steady state constraints
    if ~rem(it-1,TsParam)
        Xs = kf.x1 + Xs;
        Xsf = kff.x1 + Xsf;
        % Equilibrium constraints
        dR = CoolProp.PropsSI('D','P',Y(it,3),'H',Y(it,4),'CO2');
        Xs(3) = Xs(3)*0.05 + d1*0.95;
        Xs(7) = Xs(7)*0.05 + dR*0.95;
        Xsf(3) = Xsf(3)*0.05 + d1f*0.95;
        Xsf(7) = Xsf(7)*0.05 + dR*0.95;
        w = 0.11;
        Xs(4) = Xs(4)*0.5 + 0.5*(CoolProp.PropsSI('T','P',Y(it,1),'H',Y(it,5),'CO2')/(w+1) + w/(w+1)*TA0);
        Xsf(4) = Xsf(4)*0.5 + 0.5*(CoolProp.PropsSI('T','P',Yf(it,1),'H',Yf(it,5),'CO2')/(w+1) + w/(w+1)*TA0);
        if constrainedEstimator
            % Extreme constraints
            Xs(4) = -constrainer(-Xs(4),-TA0,2);
            Xsf(4) = -constrainer(-Xsf(4),-TA0,2);
            Xs(7) = -constrainer(-Xs(7),-80,20);
            Xsf(7) = -constrainer(-Xsf(7),-80,20);
            Xs(6) = constrainer(Xs(6),hBP+50e3,10e3);
            Xsf(6) = constrainer(Xsf(6),hBPf+50e3,10e3);
            Xs(6) = -constrainer(-Xs(6),-200e3,10e3);
            Xsf(6) = -constrainer(-Xsf(6),-200e3,10e3);
        end
    end
    % ------------ Parameter estimation -----------
    % This part is already done in the simulink file, hence clear up here.
    % The code can be found in estimatortest.m as well.
    
    % Parameter identification gets stuck in local minimum, but this would
    % be the way:
    %     DTf = TBPf - TA0;
    %     Wf = [DQf-phi(2)*rls.t(2)'; rls.t(2)]/DTf;
    % ------------ Residual Whitening and Fault Detection -----------
    % Finally whitening is not applied for simulation data
    if it > start+10
        ew = resrecord(:,it);
    else
        ew = 0;
    end
    % The EM has been already evaluated in the simulink file
    [~,fault1] = fdCUSUM.detect(ew(1));
    [~,fault2] = fdGLR.detect(ew(1));
    % ------------ Equilibrium constraints for inputs -----------
    hBP = Y(it+1,2);
    try
        TBP = CoolProp.PropsSI('T','P',Y(it+1,1),'H',hBP,'CO2');
        dBP = CoolProp.PropsSI('D','P',Y(it+1,1),'H',Y(it+1,2),'CO2');
        d1 = CoolProp.PropsSI('D','P',Y(it+1,1),'H',Y(it+1,5),'CO2');
    catch
        TBP = CoolProp.PropsSI('T','P',Y(it+1,1)+1e4,'H',hBP,'CO2');
        dBP = CoolProp.PropsSI('D','P',Y(it+1,1)+1e4,'H',Y(it+1,2),'CO2');
        d1 = CoolProp.PropsSI('D','P',Y(it+1,1)+1e4,'H',Y(it+1,5),'CO2');
    end
    U(it+1,6) = dBP;
    % In case of a fault
    if it >= detectiontime
        hBPf = outcorrecord(2,it);
        hHR = Uf(it-delayh,10);
        dG = Uf(it-delayh,7);
        CRIT = Uf(it-delayh,3);
        DmQ = Uf(it-delayh,5);
        DmV = dG*1e-4*CRIT*48 + DmQ;
        DQcor = outcorrecord(1,it);
        hBPfDQ = hHR - DQcor/DmV;
        outcorrecord(3,it) = hBPfDQ;
        try
            dBPf = CoolProp.PropsSI('D','P',Yf(it+1,1),'H',hBPf,'CO2');
        catch
            dBPf = CoolProp.PropsSI('D','P',Yf(it+1,1)+1e4,'H',hBPf,'CO2');
        end
        Yf(it+1,2) = hBPf;
        Yf(it+1,5) = Y(it+1,5) - (1-hRatio)*(hBP - hBPf);
        d1f = CoolProp.PropsSI('D','P',Yf(it+1,1),'H',Yf(it+1,5),'CO2');
        Uf(it+1,6) = dBPf;
    else
        TBPf = TBP;
        hBPf = hBP;
        Yf(it+1,2) = hBPf;
        Yf(it+1,5) = Y(it+1,5) - (1-hRatio)*(hBP - hBPf);
        d1f = d1;
        dBPf = dBP;
        Uf(it+1,6) = dBPf;
    end
    % ------------ Parameter substitutions for new iteration -----------
    if ~rem(it-1,TsParam)
        paramSampleTime = it;
        UXYW = [U(it,:)'; Xs; Y(it,1:ny)'; W];
        UXYWf = [Uf(it,:)'; Xsf; Yf(it,1:ny)'; Wf];
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
    end
    % ------------ New setpoint for linearization -----------
    if ~rem(it-1,TsParam)
        kf.x1 = zeros(nx,1);
        kff.x1 = zeros(nx,1);
        kf.xf = zeros(nx,1);
        kff.xf = zeros(nx,1);
        inputSample = U(it,:)';
        inputSamplef = Uf(it,:)';
        measSample = g(YFunction,Xs,U(it,:));
        measSamplef = g(YFunction,Xsf,Uf(it,:));
    end
    measOld = meas;
    measfOld = measf;
    meas = Y(it+1,1:ny)';
    measf = Yf(it+1,1:ny)';
    u = U(it+1,:)' - inputSample;
    uf = Uf(it+1,:)' - inputSamplef;
    y = meas - measSample;
    yf = measf - measSamplef;
    uy = [u; y; A; B; C; D; Q];
    uyf = [uf; yf; Af; Bf; Cf; Df; Qf];
    % ------------ Recording -----------
    statecorrection = kf.Kx*kf.e;
    statecorrectionf = kff.Kx*kff.e;
    traceP1 = trace(kf.P1);
    traceP1f = trace(kff.P1);
    record(:,it-start+1) = [Xs+kf.xf; traceP1; statecorrection; kf.e; W];
    recordf(:,it-start+1) = [Xsf+kff.xf; traceP1f; statecorrectionf; kff.e; Wf];
    grecord(1:2,it+1) = [fdCUSUM.g; fdGLR.g];
    faultrecord(1:2,it+1) = [fault1; fault2];
end
toc
% plotting_simulink