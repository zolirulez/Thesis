close all
addpath('C:\Users\u375749\Documents\Thesis\Codes\ModelicaPlant\bestsimulations')
addpath('C:\Users\u375749\Documents\Thesis\Codes\ModelicaPlant\objects_inits')

% Fault plot
h = figure(1);
load faultcontrol
plot(fault_sim.signals.values)
load faultcontrol_ta0
hold on
plot(fault_sim.signals.values)
hold off
xlabel('Time [s]')
ylabel('Fault [K]')
legend('Flat fault','Insulation fault')
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

% -----------------------------------------------------------------
simulationstring = 'faultcontrol_ta0';
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

% State observations
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
    '$\Delta h$, fault considered'},'Location','southwestoutside'...
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
%main_meas