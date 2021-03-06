% Script for testing the estimator on simulation or field data
% Zoltan Mark Pinter, Master Thesis, 2019

if ~exist('fielddata')
    clearvars
else
    load fielddata
end
format long
% Initializing FMIKit and adding paths
addpath('C:\Users\u375749\Documents\Thesis\Codes\Linearization')
addpath('C:\Users\u375749\Documents\Thesis\Codes\MatlabPlant')
addpath('C:\Users\u375749\Documents\Thesis\Codes\ModelicaPlant\objects_inits')
addpath('C:\Users\u375749\Documents\Thesis\Codes\ModelicaPlant\matfiles')
% 'Running modelcreation_simplified.m
modelcreation_simplified
% Running substitution_simplified.m
if ~exist('fielddata')
    [A,B,C,D] = constantsave(A,B,C,D,Qcont);
else
    addpath('C:\Users\u375749\Documents\Thesis\Codes\Measurements')
    [A,B,C,D] = fieldconstantsave(A,B,C,D,Qcont);
end
substitution_simplified


% Output function
YFunction = matlabFunction(ginit(y));
UFunction = u;
XFunctionSym = x;

% Time parameters, data sets

delayT = 300;
delayh = 50;
if ~exist('fielddata')
    % This case is outdate and is replaced by estimatortest_simulink.m
    load uy_sim_faultcontrol13 
    U = uy_sim.signals.values(:,1:nu);
    Y = uy_sim.signals.values(:,nu+1:nu+ny);
    start = 2001;
    finish = length(Y)-1;
    U(U(:,5)>0.5,5) = 0.5;
else
    start = 1001;
    finish = length(Y)-1;
    Y(:,4) = Y(:,4) - 170e3;
end
hRatio = 1/3;
Y(:,5) = hRatio*[U(round(delayh/2),10)*ones(delayh,1); U(1:end-delayh,10)] + (1-hRatio)*Y(:,2);
ny = size(Y,2);

% Running initialization functions for tools
kfinit_simulink
rlsinit_simulink
uy = [zeros(nu+ny,1); reshape(system.A,nx*nx,1); reshape(system.B,nx*nu,1); reshape(system.C,ny*nx,1); reshape(system.D,ny*nu,1); reshape(noise.Q,nx*nx,1)];
if exist('fielddata')
    Xs(1) = 60e5;
    Xs(2) = 400e3;
    Xs(3) = 250;
    Xs(4) = 40+273;
    Xs(5) = 35e5;
    Xs(6) = 300e3;
    Xs(7) = 250;
    Xs(8) = 150e3;
else
    Xs(4) = 40+273;
end

% Initial numbers
TBP = CoolProp.PropsSI('T','P',Y(start,1),'H',Y(start,2),'CO2');
dR = CoolProp.PropsSI('D','P',Y(start,3),'H',Y(start,4),'CO2');
d1 = CoolProp.PropsSI('D','P',Y(start,1),'H',Y(start,5),'CO2');
TA0 = U(start,12);
hBP = Y(start,2);
% Parameter estimation (not used now)
nw = 2;
W = [5000; 6000];

% Fault detector and operation initialization
faultinit_simulink
faultopinit

% Recording
record = NaN(nx+nx+ny+nw+1,finish-start+1);
resrecord = NaN(size(rlsInitial.t,2),finish-start+1);
paramrecord = NaN(size(rls.t,1)*size(rls.t,2),finish-start+1);
outrecord = NaN(size(rlsInitial.t,2),finish-start+1);
outcorrecord = NaN(size(rlsInitial.t,2)+1,finish-start+1);
phirecord = NaN(size(rlsInitial.t,1),finish-start+1);
grecord = NaN(3,finish-start+1);
faultrecord = NaN(3,finish-start+1);
recordf = record;

% Setpoint sampling time
TsParam = 1000;

% Further variables
meas = Y(start,:);
measf = Yf(start,:);
inputSample = U(start,:)';
inputSamplef = Uf(start,:)';
constrainedEstimator = 1;

% Noise control. Finally not applied
rng(12345)

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
    TA0 = U(it-delayT,12);
    THR = CoolProp.PropsSI('T','H',U(it-delayh,10),'P',Y(it-delayh,1),'CO2');
    if ~exist('fielddata')
        CRA = U(it-delayT,1);
    else
        CRA = 1; % It was this value all along the normal operation range..
    end
    if it > 3000 && ~faultOperation
        phiold = phirecord(:,it-start+1-250);
        outold = outrecord(:,it-start+1-250)';
        rls.regression(phiold,outold);
        phiold = phirecord(:,it-start+1-500);
        outold = outrecord(:,it-start+1-500)';
        rls.regression(phiold,outold);
    end
    CRV = U(it,2);
    CRIT = U(it,3);
    DmQ = U(it-delayh,5);
    hHR = U(it-delayh,10);
    pGC = Y(it,1);
    pR = Y(it,3);
    dG = U(it-delayh,7);
    DmIT = dG*VValues(3)*CRIT*fValues(1);
    DmV = DmQ + DmIT;
    DQ = DmV*(U(it-delayh,10) - Y(it,2));
    phi = [1; hHR; (THR-TA0); CRA; DmIT*hHR; DmQ*hHR; CRV*sqrt(pGC-pR)];
    out = [DQ hBP]; 
    rls.regression(phi,out);
    % Parameter identification gets stuck in local minimum, but this would
    % be the way:
    %     DTf = TBPf - TA0;
    %     Wf = [DQf-phi(2)*rls.t(2)'; rls.t(2)]/DTf;
    % ------------ Fault Detection -----------
    if it == start
        rls.e = zeros(1,length(rls.e));
    end
    if it > start+10
        if exist('fielddata')
            Apol = [1 -0.9945]; Bpol = [5.46]; Cpol = [1 -0.3745 -0.1045 0.2085];
            ew = ((filter(Apol,Cpol,resrecord(1,1:it-start)) - filter(Bpol,Cpol,U(start+1:it,12)-mean(U(start+1:it,12)))')')*sum(Cpol)/sum(Apol);
        else
            Apol = [1 -0.991 0.9945]; Bpol = [0.2985]; Cpol = [1 -1.852 0.8701];
            ew = (filter(Apol,Cpol,resrecord(1,1:it-start)) - filter(Bpol,Cpol,U(start+1:it,12)-mean(U(start+1:it,12)))')';
        end
    else
        ew = 0;
    end
    [~,fault3] = fdEM.detect(rls.e');
    [~,fault1] = fdCUSUM.detect(ew(end,1));
    [~,fault2] = fdGLR.detect(ew(end,1));
    % ------------ Fault Operation -----------
    if it > start + round(finish/4)
        if all(faultrecord(3,it-19-start:it-start)) && ~faultOperation
            disp('Fault detected!')
            if ~detectiontime
                detectiontime = it;
            end 
            Wsave = reshape(paramrecord(:,detectiontime-500-start),length(phi),length(out));
            faultOperation = 1;
        end
        if ~any(faultrecord(3,it-19-start:it-start))
            switchofftime = it;
        end
    end
    TauFault = TsParam;
    if detectiontime > 0 
        outcor = phi'*Wsave;
        rls.t = Wsave;
    else
        outcor = NaN;
    end
    % Equilibrium constraint values influenced by a fault
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
    if detectiontime > 0
        hBPf = outcor(2);
        hHR = Uf(it-delayh,10);
        dG = Uf(it-delayh,7);
        CRIT = Uf(it-delayh,3);
        DmQ = Uf(it-delayh,5);
        DmV = dG*1e-4*CRIT*48 + DmQ;
        DQcor = outcor(1);
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
    resrecord(:,it-start+1) = rls.e;
    paramrecord(:,it-start+1) = reshape(rls.t,size(rls.t,1)*size(rls.t,2),1);
    outrecord(:,it-start+1) = out';
    outcorrecord(1:2,it+1) = outcor';
    phirecord(:,it-start+1) = phi;
    grecord(:,it-start+1) = [fdCUSUM.g; fdGLR.g; fdEM.g];
    faultrecord(:,it-start+1) = [fault1; fault2; fault3];
end
toc