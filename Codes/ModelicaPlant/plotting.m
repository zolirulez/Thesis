t = start:Ts:finish;
figure(1)
clf
subplot(221)
hold on
plot(t,record(1,:)/1e5)
plot(t,record(5,:)/1e5)
hold off
ylim([30 100])
xlabel('Time [s]')
ylabel('Pressure [bar]')
subplot(222)
hold on
plot(t,record(2,:)/1000)
plot(t,record(6,:)/1000)
plot(t,record(8,:)/1000)
hold off
ylim([-50 600])
xlabel('Time [s]')
ylabel('Enthalpy [kJ/kg]')
subplot(223)
hold on
plot(t,record(4,:)-273.15)
hold off
ylim([-10 120])
xlabel('Time [s]')
ylabel('Temperature [C]')
subplot(224)
hold on
plot(t,record(3,:))
plot(t,record(7,:))
hold off
ylim([0 900])
xlabel('Time [s]')
ylabel('Density [kg/m^3]')

figure(2)
subplot(311)
plot(t,record(nx+1:nx+nx,:))
xlabel('Time [s]')
ylabel('Eigenvalues of P_1')
subplot(312)
plot(t,record(nx+nx+1:nx+nx+nx,:))
xlabel('Time [s]')
ylabel('State correction of K_xe')
legend('p1','h1','d1','TA1','pR','hR','dR','delta_h');
subplot(3,3,7)
plot(t,record(nx+nx+nx+1,:)/1e5,t,record(nx+nx+nx+3,:)/1e5)
legend('p_1','p_R')
xlabel('Time [s]')
ylabel('Innovation [bar]')
subplot(3,3,8)
plot(t,record(nx+nx+nx+2,:)/1e3,t,record(nx+nx+nx+4,:)/1e3)
legend('h_B_P','h_R')
xlabel('Time [s]')
ylabel('Innovation [kJ/kg]')
subplot(3,3,9)
plot(t,record(nx+nx+nx+6,:))
legend('T_A_1')
xlabel('Time [s]')
ylabel('Innovation [K]')

figure(3)
plot(t,record(nx+nx+nx+ny+1:nx+nx+nx+ny+2,:))
legend('s_0','k')
xlabel('Time [s]')
ylabel('Parameters')

figure(10)
subplot(321)
plot(start:finish,(paramrecord./paramrecord(:,end))')
grid on
xlabel('Time [s]')
ylabel('Normalized parameters')
ylim([-5 5])
subplot(312)
plot(start:finish,(resrecord./max(resrecord(:,5:end)')')')
grid on
ylim([-5 1])
xlabel('Time [s]')
ylabel('Normalized weighted residuals')
legend('DQ','DT','delta_h - used now')
subplot(313)
plot(detectiontime:length(outcorrecord)+detectiontime-1,diag([1 5e3 1])*outcorrecord',...
    start:finish,outrecord'*diag([1 5e3 1]),...
    [detectiontime detectiontime],[0 1.5e5],'--',...
    [detectiontime-300 detectiontime-300],[0 1.5e5],'--');
legend('Reestimated DQ','Reestimated DT*5e3','Reestimated delta_h',...
    'Faulty DQ','Faulty DT*5e3','Faulty delta_h','Estimated detection time','Parameter sampling')
xlabel('Time [s]')
ylabel('Outputs')
subplot(322)
fault = fault_sim.signals.values;
plot(detectiontime:length(outcorrecord)+detectiontime-1,outrecord(2,detectiontime-start:it-start)-outcorrecord(:,2)',...
    4000:finish,fault(4000:end-1)')
xlabel('Time [s]')
ylabel('Fault estimation [K]')
legend('Estimated fault','True fault')
%% Residual
global var_resid FalseAlarmTime
resid = resrecord(end,:);
resid(1) = 0;
figure(5)
resid_normal = [0 resid(1,2:2000)];
subplot(221)
autocorr(resid_normal)
subplot(222)
[pxx,f,pxxc] = periodogram(resid_normal,rectwin(length(resid_normal)),...
    length(resid_normal),1,'ConfidenceLevel',0.99);
cpxx = cumsum(pxx)./sum(pxx);
plot(f,cpxx,'b',f,f./max(f)+0.05,'r--',f+0.05*max(f),f./max(f),'r--')
xlabel('Frequency')
title('Cumulative periodrogram with 99% conf. interval')
subplot(223)
probplot('normal',resid_normal)
subplot(224)
var_resid1 = var(resid_normal);
histfit(resid_normal,50)
title(['Histogram of data with variance ' num2str(var_resid1)])
figure(6)
wcn = 0.2;
[B,A] = butter(1,wcn,'high');
resid2 = filter(B,A,resid);
resid_normal2 = resid2(1,5:2000);
subplot(221)
autocorr(resid_normal2)
subplot(222)
[pxx,f,pxxc] = periodogram(resid_normal2,rectwin(length(resid_normal2)),...
    length(resid_normal2),1,'ConfidenceLevel',0.99);
cpxx = cumsum(pxx)./sum(pxx);
plot(f,cpxx,'b',f,f./max(f)+0.05,'r--',f+0.05*max(f),f./max(f),'r--')
xlabel('Frequency')
title('Cumulative periodrogram with 99% conf. interval')
subplot(223)
probplot('normal',resid_normal2)
subplot(224)
var_resid2 = var(resid_normal2);
histfit(resid_normal2,50)
title(['Histogram of data with variance ' num2str(var_resid2)])
% Preparation for detection
[B,A] = butter(3,0.9*wcn,'low');
resid = filter(B,A,resid);
resid_normal = [0 resid(1,2:2000)];
var_resid = var(resid_normal);
% --------------- Detection ----------------------
figure(7)
subplot(311)
plot(start:finish,resid')
grid on
ylabel('Residual')
xlabel('Time [s]')
ylim([-2e4 2e4])
subplot(334)
plot(start:finish,grecord(1,:)',start:finish,fdCUSUM.h*rectwin(length(resid)),'r--')
ylabel('g function')
xlabel('Time [s]')
title('CUSUM')
subplot(337)
plot(start:finish,grecord(1,:)'>fdCUSUM.h,'LineWidth',2)
ylabel('Fault detection')
xlabel('Time [s]')
% GLR
figure(7)
subplot(335)
plot(start:finish,grecord(2,:)',start:finish,fdGLR.h*rectwin(length(resid)),'r--')
ylim([0 fdGLR.h*10])
ylabel('g function')
xlabel('Time [s]')
title('GLR')
subplot(338)
plot(start:finish,grecord(2,:)'>fdGLR.h,'LineWidth',2)
ylabel('Fault detection')
xlabel('Time [s]')
% Moving square change detection
figure(7)
M = 500;
subplot(336)
plot(start:finish,grecord(3,:)',start:finish,fdEM.h*rectwin(length(resid)),'r--')
ylabel('g function')
xlabel('Time [s]')
title('Expectation Maximization (multidim.)')
subplot(339)
plot(start:finish,grecord(3,:)'>fdEM.h,'LineWidth',2)
ylabel('Fault detection')
xlabel('Time [s]')
