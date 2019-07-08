% Script for testing the estimator on simulation or field data

if ~exist('fielddata')
    clearvars
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
DXFunction = matlabFunction(finit(Dx));
UFunction = u;
XFunctionSym = x;
load uy_sim_faulty_chirp2 % uy_sim_faulty

% Time parameters, data sets

delayT = 1*3;
delayh = 25*3;
if ~exist('fielddata')
    U = uy_sim.signals.values(:,1:nu);
    Y = uy_sim.signals.values(:,nu+1:nu+ny);
    start = 2001;
    finish = 10000;
else
    start = 501;
    finish = length(Y)-1;
end
Y(:,5) = 1/3*[U(round(delayh/2),10)*ones(delayh,1); U(1:end-delayh,10)] + 2/3*Y(:,2);
ny = size(Y,2);

% Running initialization functions for tools
kfinit_simulink
rlsinit_simulink
uy = [zeros(nu+ny,1); reshape(system.A,nx*nx,1); reshape(system.B,nx*nu,1); reshape(system.C,ny*nx,1); reshape(system.D,ny*nu,1); reshape(noise.Q,nx*nx,1)];
% if exist('fielddata')
%     Xs(1) = 60e5;
%     Xs(2) = 400e3;
%     Xs(3) = 250;
%     Xs(4) = 60+273;
%     Xs(5) = 35e5;
%     Xs(6) = 300e3;
%     Xs(7) = 250;
%     Xs(8) = 150e3;
%     Xs(9) = 0;
% else
Xs(4) = 60+273;
Xsf(4) = Xs(4);
% end

% Initial numbers
TBP = CoolProp.PropsSI('T','P',Y(start,1),'H',Y(start,2),'CO2');
dR = CoolProp.PropsSI('D','P',Y(start,3),'H',Y(start,4),'CO2');
d1 = CoolProp.PropsSI('D','P',Y(start,1),'H',Y(start,5),'CO2');
TA0 = U(start,12);
% Parameter estimation (not used now)
nw = 2;
W = [5000; 6000];

% Fault detector and operation initialization
faultinit_simulink

% Recording
record = NaN(nx+nx+nx+ny+nw,finish-start+1);
resrecord = NaN(size(rlsInitial.t,2),finish-start+1);
paramrecord = NaN(size(rls.t,1)*size(rls.t,2),finish-start+1);
outrecord = NaN(size(rlsInitial.t,2),finish-start+1);
outcorrecord = NaN(size(rlsInitial.t,2),finish-start+1);
phirecord = NaN(size(rlsInitial.t,1),finish-start+1);
grecord = NaN(3,finish-start+1);
faultrecord = NaN(3,finish-start+1);
recordf = record;

% For earlier states
% DQo = 74e3;
% DQ = 74e3;
% DTo = 3;
% DT = 3;
meas = Y(start,:);
measf = Yf(start,:);

for it = start:finish
    if ~rem(it,100)
        disp(['Iteration ' num2str(it)])
    end
    % ------------ State estimation -----------
    kf.markovPredictor(uy(1:nu),uy(nu+1:nu+ny),reshape(uy(nu+ny+1:nu+ny+nx*nx),nx,nx),reshape(uy(nu+ny+nx*nx+1:nu+ny+nx*nx+nx*nu),nx,nu),reshape(uy(nu+ny+nx*nx+nx*nu+1:nu+ny+nx*nx+nx*nu+ny*nx),ny,nx), reshape(uy(nu+ny+nx*nx+nx*nu+ny*nx+1:nu+ny+nx*nx+nx*nu+ny*nx+ny*nu),ny,nu), reshape(uy(nu+ny+nx*nx+nx*nu+ny*nx+ny*nu+1:nu+ny+nx*nx+nx*nu+ny*nx+ny*nu+nx*nx),nx,nx));
    kff.markovPredictor(uyf(1:nu),uyf(nu+1:nu+ny),reshape(uyf(nu+ny+1:nu+ny+nx*nx),nx,nx),reshape(uyf(nu+ny+nx*nx+1:nu+ny+nx*nx+nx*nu),nx,nu),reshape(uyf(nu+ny+nx*nx+nx*nu+1:nu+ny+nx*nx+nx*nu+ny*nx),ny,nx), reshape(uyf(nu+ny+nx*nx+nx*nu+ny*nx+1:nu+ny+nx*nx+nx*nu+ny*nx+ny*nu),ny,nu), reshape(uyf(nu+ny+nx*nx+nx*nu+ny*nx+ny*nu+1:nu+ny+nx*nx+nx*nu+ny*nx+ny*nu+nx*nx),nx,nx));
    %     Xs = Xs + Ts*f(DXFunction,Xs,U(it,:));
    %     Xsf = Xsf + Ts*f(DXFunction,Xsf,U(it,:));
    % Steady state constraints on unobservable density
    if ~rem(it-1,100)
        Xsnew = kf.x1 + Xs;
        Xsfnew = kff.x1 + Xsf;
        Xsnew(3) = Xsnew(3)*0.4 + d1*0.6;
        Xsnew(7) = Xsnew(7)*0.4 + dR*0.6;
        Xsfnew(3) = Xsfnew(3)*0.4 + d1f*0.6;
        Xsfnew(7) = Xsfnew(7)*0.4 + dR*0.6;
        w = 0.1;
        Xsnew(4) = Xsnew(4)*0.4 + 0.6*(CoolProp.PropsSI('T','P',Xsnew(1),'H',Xsnew(2),'CO2')/(w+1) + w/(w+1)*TA0);
        Xsfnew(4) = Xsfnew(4)*0.4 + 0.6*(CoolProp.PropsSI('T','P',Xsfnew(1),'H',Xsfnew(2),'CO2')/(w+1) + w/(w+1)*TA0);
    end
    Xs = 0.99*Xs + 0.01*Xsnew;
    Xsf = 0.99*Xsf + 0.01*Xsfnew;
    % ------------ Parameter estimation -----------
    TA0 = U(it-delayT,12);
    THR = CoolProp.PropsSI('T','H',U(it-delayh,10),'P',Y(it-delayh,1),'CO2');
    CRG = U(it-delayh,4);
    CRIT = U(it-delayh,3);
    CRA = U(it-delayT,1);
    DmV = U(it,2)*KvValues(1)*sqrt(U(it,6)*(Y(it,1) - Y(it,3)));
    DQ = DmV*(U(it-delayh,10) - Y(it,2));
    CRV = U(it,2);
    phi = [1; CRA; CRV; TA0-273; CRIT; CRG; THR-273];
    out = [DQ TBP]; 
    rls.regression(phi,out);
    % Parameter identification gets stuck in local minimum, but this would
    % be the way:
    %     DTf = TBPf - TA0;
    %     Wf = [DQf-phi(2)*rls.t(2)'; rls.t(2)]/DTf;
    % ------------ Fault Detection -----------
    if it == start
        rls.e = zeros(1,length(rls.e));
    end
    [~,fault3] = fdEM.EM(rls.e');
    if it > start+10
        if exist('fielddata')
%             Apol = [1 -0.9886]; Bpol = [-1.078]; Cpol = [1 -0.3908 0.04806 0.3148];
            Apol = [1 -0.991]; Bpol = [-0.8966]; Cpol = [1 -0.424 0.05303 0.2916];
            ew = filter(Apol,Cpol,resrecord(1,1:it-start)) - filter(Bpol,Cpol,U(start+1:it,12)-mean(U(start+1:it,12)))';
        else
            ew = rls.e; %filter(Apol,Cpol,resrecord(1,1:it-start)) ;
        end
    else
        ew = 0;
    end
    [~,fault1] = fdCUSUM.CUSUM(ew(1,end));
    [~,fault2] = fdGLR.GLR(ew(1,end));
    % ------------ Fault Operation -----------
    if it > start + 20
        if all(faultrecord(3,it-19-start:it-start)) && ~faultOperation
            Wsave = reshape(paramrecord(:,it-500-start),length(phi),length(out));
            if ~detectiontime
                detectiontime = it;
            end
            faultOperation = 1;
        end
        if ~any(faultrecord(3,it-19-start:it-start))
            faultOperation = 0;
            switchofftime = it;
        end
    end
    TauFault = param.ResponsibilityTimeConstant*10;
    if faultOperation
        outcor = phi'*Wsave;
        Tfault = (1-Ts/TauFault)*Tfault + Ts/TauFault*(out(2) - outcor(2));
        rls.t = Wsave;
    else
        Tfault = 0.5*Tfault;
        outcor = NaN;
    end
    % ------------ Parameter substitutions for new iteration -----------
    if ~rem(it-1,100)
        paramSampleTime = it;
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
    end
    % ------------ Constraints -----------
    hBP = Y(it+1,2);
    TBP = CoolProp.PropsSI('T','P',Y(it+1,1),'H',hBP,'CO2');
    dR = CoolProp.PropsSI('D','P',Y(it+1,3),'H',Y(it+1,4),'CO2');
    d1 = CoolProp.PropsSI('D','P',Y(it+1,1),'H',Y(it+1,5),'CO2');
    % In case of a fault
    TBPf = TBP - Tfault;
    hBPf = CoolProp.PropsSI('H','P',Yf(it+1,1),'T',TBPf,'CO2');
    Yf(it+1,2) = hBPf;
    Yf(it+1,5) = Y(it+1,5) - 2/3*(hBP - hBPf);
    d1f = CoolProp.PropsSI('D','P',Yf(it+1,1),'H',Yf(it+1,5),'CO2');
    dBPf = CoolProp.PropsSI('D','P',Yf(it+1,1),'H',Yf(it+1,2),'CO2');
    Uf(it+1,6) = dBPf;
    % ------------ New setpoint for linearization -----------
    if ~rem(it-1,100)
        kf.x1 = zeros(nx,1);
        kff.x1 = zeros(nx,1);
        inputSample = U(it,:)';
        inputSamplef = Uf(it,:)';
        measSample = g(YFunction,Xsnew,U(it,:));
        measSamplef = g(YFunction,Xsfnew,Uf(it,:));
    end
    measOld = meas;
    measfOld = measf;
    meas = Y(it+1,1:ny)';
    measf = Yf(it+1,1:ny)';
%     meas(2) = 0.01*meas(2) + 0.99*measOld(2);
%     meas(5) = 0.01*meas(5) + 0.99*measOld(5);
%     measf(2) = 0.01*measf(2) + 0.99*measfOld(2);
%     measf(5) = 0.01*measf(5) + 0.99*measfOld(5);
    u = U(it+1,:)' - inputSample;
    uf = Uf(it+1,:)' - inputSamplef;
    y = meas - measSample;
    yf = measf - measSamplef;
    uy = [u; y; A; B; C; D; Q];
    uyf = [uf; yf; Af; Bf; Cf; Df; Qf];
    % ------------ Recording -----------
    statecorrection = kf.Kx*kf.e;
    statecorrectionf = kff.Kx*kff.e;
    eigP1 = eig(kf.P1); 
    eigP1f = eig(kff.P1); 
    record(:,it-start+1) = [Xs+kf.xf; eigP1; statecorrection; kf.e; W];
    recordf(:,it-start+1) = [Xsf+kff.xf; eigP1f; statecorrectionf; kff.e; Wf];
    resrecord(:,it-start+1) = rls.e;
    paramrecord(:,it-start+1) = reshape(rls.t,size(rls.t,1)*size(rls.t,2),1);
    outrecord(:,it-start+1) = out';
    outcorrecord(:,it-start+1) = outcor';
    phirecord(:,it-start+1) = phi;
    grecord(:,it-start+1) = [fdCUSUM.g; fdGLR.g; fdEM.g];
    faultrecord(:,it-start+1) = [fault1; fault2; fault3];
end

plotting