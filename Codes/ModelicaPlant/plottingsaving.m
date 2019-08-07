clearvars
close all
addpath('C:\Users\u375749\Documents\Thesis\Codes\ModelicaPlant\bestsimulations')
addpath('C:\Users\u375749\Documents\Thesis\Codes\ModelicaPlant\objects_inits')

% Fault plot
h = figure(1);
load faultcontrol
faultconversion
plot(faultest(:,1))
hold on
plot(fault_sim.signals.values,'LineWidth',1.5)
load faultcontrol_ta0
faultconversion
plot(faultest(:,1))
plot(fault_sim.signals.values,'LineWidth',1.5)
hold off
xlabel('Time [s]')
ylabel('Fault [K]')
legend({'Estimation of:','\quad flat fault','Estimation of:',...
    '\quad insulation fault'},'location','southeast','Interpreter','latex')
saveas(h,'fault.png')

% COP plot
h = figure(2);
load faultignore
plot(COP_sim.signals.values(:,1))
load faultignore_ta0
hold on
plot(COP_sim.signals.values(:,1))
load faultcontrol
plot(COP_sim.signals.values(:,1),'LineWidth',2)
load faultcontrol_ta0
plot(COP_sim.signals.values(:,1))
hold off
xlabel('Time [s]')
ylabel('Coefficient of Performance [-]')
xlim([1000 20000])
legend('Ignored flat fault','Ignored insulation fault',...
    'Fault operation, flat fault','Fault operation, insulation fault')
saveas(h,'COP.png')

% Input saturation plot
h = figure(3);
load faultignore_ta0
plot(inputs_sim.Data(:,[1 5]))
load faultcontrol_ta0
hold on
plot(inputs_sim.Data(:,[1 5]))
hold off
xlabel('Time [s]')
ylabel('Input capacity ratios [-]')
xlim([1000 20000])
ylim([0 1.1])
legend({'$CR_V$, ignored fault,','$CR_V$, fault operation',...
    '$CR_{IT}$, ignored fault,','$CR_{IT}$, fault operation'},...
    'Interpreter','latex')
saveas(h,'inputs.png')

% Ambient temperature plot
h = figure(4);
plot(uy_sim.signals.values(:,12)-273.15)
xlabel('Time [s]')
ylabel('Ambient temperature [C]')
xlim([1000 20000])
saveas(h,'ambient.png')

% Enthalpy GC and evaporator plot
h = figure(5);
load faultignore_ta0
plot(446.5e3*rectwin(20000),'LineWidth',2)
hold on
plot(Dm_h_sim.signals.values(:,[6 7])/1e3)
load faultcontrol_ta0
plot(Dm_h_sim.signals.values(:,[6 7])/1e3,'--')
hold off
xlabel('Time [s]')
ylabel('Evaporator outlet enthalpies [kJ/kg]')
xlim([1000 20000])
ylim([430 470])
legend({'Reference value','Cooler, ignored fault,','Freezer, ignored fault',...
    'Cooler, fault operation,','Freezer, fault operation'})
saveas(h,'evapenthalpies.png')


% GC entahlpy plot
h = figure(6);
load TGC_faultignore_ta0
plot(hcw/1e3,'LineWidth',3)
hold on
load faultignore_ta0
plot(h1_sim.signals.values(:,end-9:end)/1e3)
hold off
xlabel('Time [s]')
ylabel('Gas cooler cell enthalpies [kJ/kg]')
xlim([1000 20000])
ylim([250 500])
legend('Weighted average')
saveas(h,'GCenthalpies.png')

% -----------------------------------------------------------------
simulationstring = 'faultignore_ta0';
TsParam = 100;
estimatortest_simulink

% Measurement reestimation
h = figure(7);
plot(1000:finish,Y(1001:end,2)/1000)
hold on
plot(1000:finish,outcorrecord(4,1001:end)/1000,'LineWidth',2)
plot(1000:finish,outcorrecord(2,1001:end)/1000)
plot(1000:finish,outcorrecord(3,1001:end)/1000)
plot([5000 5000],[0 1.5e5],'--')
hold off
ylabel('h_B_P [kJ/kg]')
xlabel('Time [s]')
legend({'$h_{BPm}$','$h_{BP}$','$\hat{h}_{BP}$',...
    '$h_{HR}-\hat{Q}/\dot{m}_{HR}$',...
    'Fault occurence'},'Location','southeast'...
    ,'Interpreter','Latex');
ylim([100 450])
xlim([1000 length(Y)])
saveas(h,'hBP_faultignore_ta0.png')

% State observations
t = start:Ts:finish;
tcw = start:Ts:it;
h = figure(12);
subplot(221)
plot(tcw,hcw(start:it)/1000,'LineWidth',2)
hold on
plot(t,record(2,:)/1000,'LineWidth',1.5)
plot(t,record(8,:)/1000,'LineWidth',1.5)
plot(t,recordf(2,:)/1000,'--')
plot(t,recordf(8,:)/1000,'--')
hold off
ylim([-50 450])
xlabel('Time [s]')
ylabel('Enthalpy [kJ/kg]')
title('Fault ignored, parameter sample time: 100 s')
% saveas(h,'hhat.png')

h = figure(13);
subplot(221)
plot(tcw,Tcw(start:it)-273.15,'LineWidth',2)
hold on
plot(t,record(4,:)-273.15,'LineWidth',1.5)
plot(t,recordf(4,:)-273.15,'--')
hold off
xlabel('Time [s]')
ylabel('Temperature [C]')
title('Fault ignored, parameter sample time: 100 s')
% saveas(h,'TAhat.png')

% -----------------------------------------------------------------
simulationstring = 'faultignore_ta0';
TsParam = 1000;
estimatortest_simulink

% State observations
h = figure(12);
subplot(222)
plot(tcw,hcw(start:it)/1000,'LineWidth',2)
hold on
plot(t,record(2,:)/1000,'LineWidth',1.5)
plot(t,record(8,:)/1000,'LineWidth',1.5)
plot(t,recordf(2,:)/1000,'--')

plot(t,recordf(8,:)/1000,'--')
hold off
ylim([-50 450])
xlabel('Time [s]')
ylabel('Enthalpy [kJ/kg]')
title('Fault ignored, parameter sample time: 1000 s')
% saveas(h,'hhat.png')

h = figure(13);
subplot(222)
plot(tcw,Tcw(start:it)-273.15,'LineWidth',2)
hold on
plot(t,record(4,:)-273.15,'LineWidth',1.5)
plot(t,recordf(4,:)-273.15)
hold off
xlabel('Time [s]')
ylabel('Temperature [C]')
title('Fault ignored, parameter sample time: 1000 s')
% saveas(h,'TAhat.png')

h = figure(15);
plot(t,record(6,:)/1000,'LineWidth',1.5)
hold on
plot(t,recordf(6,:)/1000,'--','LineWidth',1.5)
hold off
xlim([1000 length(Y)])
xlabel('Time [s]')
ylabel('Receiver enthalpy [C]')
legend('Fault ignored','Fault considered')
saveas(h,'hR.png')

% Residual statistics plot
h = figure(16);
resid = resrecord(1,2001:5000)';
subplot(221)
autocorr(detrend(resid))
subplot(222)
[pxx,f,pxxc] = periodogram(resid,rectwin(length(resid)),...
    length(resid),1,'ConfidenceLevel',0.99);
cpxx = cumsum(pxx)./sum(pxx);
plot(f,cpxx,'b',f,f./max(f)+0.05,'r--',f+0.05*max(f),f./max(f),'r--')
xlabel('Frequency')
title('Cum. periodrogram with 99% conf. interval')
subplot(223)
probplot('normal',resid)
title('Prob. plot for normal dist.')
subplot(224)
var_resid1 = var(resid);
histfit(resid,50)
ylabel('Sample frequency')
title('Histogram of data')
saveas(h,'resstat.png')

% -----------------------------------------------------------------
simulationstring = 'faultcontrol_ta0';
load TGC_faultcontrol_ta0.mat
TsParam = 100;
estimatortest_simulink

% Measurement reestimation
h = figure(8);
plot(1000:finish,Y(1001:end,2)/1000)
hold on
plot(1000:finish,outcorrecord(4,1001:end)/1000,'LineWidth',2)
plot(1000:finish,outcorrecord(2,1001:end)/1000)
plot(1000:finish,outcorrecord(3,1001:end)/1000)
plot([detectiontime detectiontime],[0 1.5e5],'--',...
    [detectiontime-500 detectiontime-500],[0 1.5e5],'--',...
    [5000 5000],[0 1.5e5],'--')
hold off
ylabel('h_B_P [kJ/kg]')
xlabel('Time [s]')
legend({'$h_{BPm}$','$h_{BP}$','$\hat{h}_{BP}$',...
    '$h_{HR}-\hat{Q}/\dot{m}_{HR}$','Detection time','Parameter sampling',...
    'Fault occurence'},'Location','southeast'...
    ,'Interpreter','Latex');
ylim([230 330])
xlim([1000 length(Y)])
saveas(h,'hBP_faultcontrol_ta0.png')

% Lissajous figure
h = figure(9);
clf
subplot(121)
hold on
plot(smooth(Tcw(start+500:it)-273.15-3,100),smooth(record(4,501:end)-273.15-3,100))
plot(smooth(Tcw(start+500:it)-273.15-3,100),smooth(recordf(4,501:end)-273.15-3,100))
plot(linspace(10,60,2),linspace(10,60,2))
hold off
xlabel('Theoretical approx. of temperature [C]')
ylabel('Estimation [C]')
legend('Original','Reestimated')
title('Parameter sampling time: 100 s')
xlim([15 55])
ylim([15 55])
hold off

% Trace of P1 and innovation
h = figure(10);
subplot(211)
plot(t,record(nx+1,:))
xlabel('Time [s]')
ylabel('Trace of P_1')
xlim([1000 length(Y)])
subplot(234)
plot(t,record(nx+nx+2,:)/1e5,t,record(nx+nx+4,:)/1e5)
legend('p_1','p_R')
xlabel('Time [s]')
ylabel('Innovation [bar]')
xlim([1000 length(Y)])
subplot(235)
plot(t,record(nx+nx+3,:)/1e3,t,record(nx+nx+6,:)/1e3)
legend('h_B_P','h_1')
xlabel('Time [s]')
ylabel('Innovation [kJ/kg]')
xlim([1000 length(Y)])
subplot(236)
plot(t,record(nx+nx+5,:)/1e3)
legend('h_R')
xlabel('Time [s]')
ylabel('Innovation [kJ/kg]')
xlim([1000 length(Y)])
saveas(h,'innotrace.png')

% State observations: enthalpy
h = figure(12);
subplot(223)
plot(tcw,hcw(start:it)/1000,'LineWidth',2)
hold on
plot(t,record(2,:)/1000,'LineWidth',1.5)
plot(t,record(8,:)/1000,'LineWidth',1.5)
plot(t,recordf(2,:)/1000,'--')
plot(t,recordf(8,:)/1000,'--')
hold off
ylim([-50 450])
xlabel('Time [s]')
ylabel('Enthalpy [kJ/kg]')
title('Fault operation, parameter sample time: 100 s')
% saveas(h,'hhat.png')

% State observations: temperature
h = figure(13);
subplot(223)
plot(tcw,Tcw(start:it)-273.15,'LineWidth',2)
hold on
plot(t,record(4,:)-273.15,'LineWidth',1.5)
plot(t,recordf(4,:)-273.15,'--')
hold off
xlabel('Time [s]')
ylabel('Temperature [C]')
title('Fault operation, parameter sample time: 100 s')
% saveas(h,'TAhat.png')

% Residuals and detectors
h= figure(17);
subplot(311)
plot(1:finish,resrecord(:,1:finish)/1e3)
ylabel('Residuals')
xlabel('Time [s]')
xlim([2000 length(Y)])
legend({'$\dot{\tilde{Q}}$ [kW]','$\tilde{h}_{BP}$ [kJ/kg]'},'Interpreter','latex')
subplot(334)
plot(2000:finish,grecord(1,2001:end)',0:finish,fdCUSUM.h*rectwin(length(Y)),'r--')
ylabel('Decision function, g')
xlabel('Time [s]')
title('CUSUM')
xlim([2000 length(Y)])
subplot(337)
plot(2000:finish,grecord(1,2001:end)'>fdCUSUM.h,'LineWidth',2)
ylabel('Fault detection')
xlabel('Time [s]')
xlim([2000 length(Y)])
% GLR
subplot(335)
plot(2000:finish,grecord(2,2001:end)',0:finish,fdGLR.h*rectwin(length(Y)),'r--')
ylim([0 fdGLR.h*10])
ylabel('Decision function, g')
xlabel('Time [s]')
title('GLR')
xlim([2000 length(Y)])
subplot(338)
plot(2000:finish,grecord(2,2001:end)'>fdGLR.h,'LineWidth',2)
ylabel('Fault detection')
xlabel('Time [s]')
xlim([2000 length(Y)])
% EM
M = 500;
subplot(336)
plot(2000:finish,grecord(3,2001:end)',0:finish,fdEM.h*rectwin(length(Y)),'r--')
ylabel('Decision function, g')
xlabel('Time [s]')
title('Expectation Maximization (multidim.)')
xlim([2000 length(Y)])
subplot(339)
plot(2000:finish,grecord(3,2001:end)'>fdEM.h,'LineWidth',2)
ylabel('Fault detection')
xlabel('Time [s]')
xlim([2000 length(Y)])
saveas(h,'resid.png')

% -----------------------------------------------------------------
simulationstring = 'faultcontrol_ta0';
TsParam = 1000;
estimatortest_simulink

% Lissajous figure
h = figure(9);
subplot(122)
plot(smooth(Tcw(start+500:it)-273.15-3,100),smooth(record(4,501:end)-273.15-3,100))
hold on
plot(smooth(Tcw(start+500:it)-273.15-3,100),smooth(recordf(4,501:end)-273.15-3,100))
plot(linspace(10,60,2),linspace(10,60,2))
hold off
xlabel('Theoretical approx. of temperature [C]')
ylabel('Estimation [C]')
legend('Original','Reestimated')
title('Parameter sampling time: 1000 s')
xlim([15 55])
ylim([15 55])
hold off
saveas(h,'lissajous_1000.png')

% State observations
h = figure(12);
subplot(224)
plot(tcw,hcw(start:it)/1000,'LineWidth',2)
hold on
plot(t,record(2,:)/1000,'LineWidth',1.5)

plot(t,record(8,:)/1000,'LineWidth',1.5)
plot(t,recordf(2,:)/1000,'--')

plot(t,recordf(8,:)/1000,'--')
hold off
ylim([-50 450])
xlabel('Time [s]')
ylabel('Enthalpy [kJ/kg]')
legend({'Theoretical approximation of $h_{GC}$','$\hat{h}_{GC}$, fault ignored',...
    '$\Delta h$, fault ignored',...
    '$\hat{h}_{GC}$, fault considered',...
    '$\Delta h$, fault considered'},'Location','southwest'...
    ,'Interpreter','Latex');
title('Fault operation, parameter sample time: 1000 s')
saveas(h,'hhat.png')

h = figure(13);
subplot(224)
plot(tcw,Tcw(start:it)-273.15,'LineWidth',2)
hold on
plot(t,record(4,:)-273.15,'LineWidth',1.5)
plot(t,recordf(4,:)-273.15,'--')
hold off
xlabel('Time [s]')
ylabel('Temperature [C]')
legend({'Theoretical approximation of $T_{GC}$','$\hat{T}_{A,GC}$, fault ignored',......
    '$\hat{T}_{A,GC}$, fault considered'},'Location','southwest'...
    ,'Interpreter','Latex');
title('Fault operation, parameter sample time: 1000 s')
saveas(h,'TAhat.png')

% -----------------------------------------------------------------
addpath('C:\Users\u375749\Documents\Thesis\Codes\Measurements')
load datatreatment
h = figure(14);
plot(hBPoriginal/1e3,'LineWidth',1.5)
hold on
plot(hBP/1e3,'--','LineWidth',1.5)
hold off
xlabel('Time [s]')
ylabel('Enthalpy measurement [kJ/kg]')
legend({'Before correction','After placing on the saturation line'},'Location','best')
saveas(h,'datatreatment.png')

% Running the field experiment
TsParam = 1000;
main_meas
faultconversion

% Fault and input saturation 7748 12367
h = figure(18);
subplot(121)
plot([zeros(1,7747) -5*ones(1,12367-7747) zeros(1,length(Y)-12368)],'LineWidth',1.5)
hold on
plot(faultest)
hold off
xlabel('Time [s]')
ylabel('Fault [K]')
legend({'Fault','Estimation from $\hat{h}_{BP}$',...
    'Estimation from $h_{HR}-\dot{\hat{Q}}/\dot{m}_{HR}$'},...
    'location','west','Interpreter','latex')
ylim([-7 2])
subplot(122)
plot(1:length(Y),U(:,[1 2 4]))
load CRMT
hold on
plot(1:length(Y), CRMT)
hold off
xlabel('Time [s]')
ylabel('Input capacity ratios [-]')
ylim([0 1.1])
legend({'$CR_A$','$CR_{V}$','$CR_{G}$','$CR_{MT}$'},...
    'location','west','Interpreter','latex')
saveas(h,'faultinputs_field.png')

% Residuals and detectors
h= figure(19);
subplot(321)
plot(start+2:length(Y),ew/1e3)
hold on
plot(start+1:length(Y),resrecord(1,:)/1e3,'LineWidth',2)
hold off
ylabel('Residuals')
xlabel('Time [s]')
xlim([2000 length(Y)])
ylim([-20 40])
legend({'Whitened $\dot{\tilde{Q}}$','$\dot{\tilde{Q}}$ [kW]'},'Interpreter','latex')
subplot(322)
plot(start+1:length(Y),resrecord(2,:)/1e3)
xlim([2000 length(Y)])
xlabel('Time [s]')
legend({'$\tilde{h}_{BP}$ [kJ/kg]'},'Interpreter','latex')
subplot(334)
plot(start:finish,grecord(1,:)',0:finish,fdCUSUM.h*rectwin(length(Y)),'r--')
ylabel('Decision function, g')
xlabel('Time [s]')
title('CUSUM')
xlim([2000 length(Y)])
subplot(337)
plot(start:finish,grecord(1,:)'>fdCUSUM.h,'LineWidth',2)
ylabel('Fault detection')
xlabel('Time [s]')
xlim([2000 length(Y)])
% GLR
subplot(335)
plot(start:finish,grecord(2,:)',0:finish,fdGLR.h*rectwin(length(Y)),'r--')
ylim([0 fdGLR.h*10])
ylabel('Decision function, g')
xlabel('Time [s]')
title('GLR')
xlim([2000 length(Y)])
subplot(338)
plot(start:finish,grecord(2,:)'>fdGLR.h,'LineWidth',2)
ylabel('Fault detection')
xlabel('Time [s]')
xlim([2000 length(Y)])
% EM
M = 500;
subplot(336)
plot(start:finish,grecord(3,:)',0:finish,fdEM.h*rectwin(length(Y)),'r--')
ylabel('Decision function, g')
xlabel('Time [s]')
title('Expectation Maximization (multidim.)')
xlim([2000 length(Y)])
subplot(339)
plot(start:finish,grecord(3,:)'>fdEM.h,'LineWidth',2)
ylabel('Fault detection')
xlabel('Time [s]')
xlim([2000 length(Y)])
saveas(h,'resid_field.png')

% State estimations
t = start:Ts:finish;
figure(20)
clf
subplot(221)
hold on
plot(t,record(1,:)/1e5,'LineWidth',1.5)
plot(t,record(5,:)/1e5,'LineWidth',1.5)
plot(t,recordf(1,:)/1e5,'--')
plot(t,recordf(5,:)/1e5,'--')
hold off
ylim([20 80])
xlabel('Time [s]')
ylabel('Pressure [bar]')
legend({'$p_{GC}$','$p_{R}$','$p_{GCc}$','$p_{Rc}$'},...
    'location','southwest','Interpreter','latex')
subplot(222)
hold on
plot(t,record(2,:)/1000,'LineWidth',1.5)
plot(t,record(6,:)/1000,'LineWidth',1.5)
plot(t,record(8,:)/1000,'LineWidth',1.5)
plot(t,recordf(2,:)/1000,'--')
plot(t,recordf(6,:)/1000,'--')
plot(t,recordf(8,:)/1000,'--')
hold off
ylim([-50 600])
xlabel('Time [s]')
ylabel('Enthalpy [kJ/kg]')
legend({'$h_{GC}$','$h_{R}$','$\Delta h$','$h_{GCc}$','$h_{Rc}$',...
    '$\Delta h_c$'},'location','southwest','Interpreter','latex')
subplot(223)
hold on
plot(t,record(4,:)-273.15,'LineWidth',1.5)
plot(t,recordf(4,:)-273.15,'--')
hold off
% ylim([-10 120])
xlabel('Time [s]')
ylabel('Temperature [C]')
legend({'$T_{A}$','$T_{Ac}$'},'location','southwest','Interpreter','latex')
subplot(224)
hold on
plot(t,record(3,:),'LineWidth',1.5)
plot(t,record(7,:),'LineWidth',1.5)
plot(t,recordf(3,:),'--')
plot(t,recordf(7,:),'--')
hold off
ylim([0 900])
xlabel('Time [s]')
ylabel('Density [kg/m^3]')
legend({'$\rho_{GC}$','$\rho_{R}$','$\rho_{GCc}$','$\rho_{Rc}$'},...
    'location','southwest','Interpreter','latex')
saveas(h,'states_field.png')

% TODO: signal whitening plots
% structure order choice plots