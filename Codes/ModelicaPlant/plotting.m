load fault_chirp
load TGC2

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
% ylim([-10 120])
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
legend('p1','h1','d1','h1','pR','hR','dR','delta_h');
subplot(3,3,7)
plot(t,record(nx+nx+nx+1,:)/1e5,t,record(nx+nx+nx+3,:)/1e5)
legend('p_1','p_R')
xlabel('Time [s]')
ylabel('Innovation [bar]')
grid on
subplot(3,3,8)
plot(t,record(nx+nx+nx+2,:)/1e3,t,record(nx+nx+nx+5,:)/1e3)
legend('h_B_P','h_1')
xlabel('Time [s]')
ylabel('Innovation [kJ/kg]')
grid on
subplot(3,3,9)
plot(t,record(nx+nx+nx+4,:)/1e3)
legend('h_R')
xlabel('Time [s]')
ylabel('Innovation [kJ/kg]')
grid on

figure(3)
plot(t,record(nx+nx+nx+ny+1:nx+nx+nx+ny+2,:))
legend('s_0','k')
xlabel('Time [s]')
ylabel('Parameters')

figure(1)
load hcw_chirp
subplot(221)
hold on
plot(t,recordf(1,:)/1e5,'--')
plot(t,recordf(5,:)/1e5,'--')
hold off
ylim([30 100])
xlabel('Time [s]')
ylabel('Pressure [bar]')
subplot(222)
hold on
plot(t,recordf(2,:)/1000,'--')
plot(t,recordf(6,:)/1000,'--')
plot(t,recordf(8,:)/1000,'--')
plot(t,hcw(start:it)/1000,'--')
hold off
ylim([-50 600])
xlabel('Time [s]')
ylabel('Enthalpy [kJ/kg]')
subplot(223)
hold on
plot(t,recordf(4,:)-273.15,'--')
plot(t,Tcw(start:it)-273)
hold off
% ylim([-10 120])
xlabel('Time [s]')
ylabel('Temperature [C]')
subplot(224)
hold on
plot(t,recordf(3,:),'--')
plot(t,recordf(7,:),'--')
hold off
ylim([0 900])
xlabel('Time [s]')
ylabel('Density [kg/m^3]')


% figure(3)
% hold on
% plot(t,recordf(nx+nx+nx+ny+1:nx+nx+nx+ny+2,:),'--')
% hold off
% legend('s_0','k')
% xlabel('Time [s]')
% ylabel('Parameters')

figure(10)
subplot(321)
plot(start:finish,(paramrecord./paramrecord(:,it-start))')
grid on
xlabel('Time [s]')
ylabel('Normalized parameters')
ylim([-20 20])
subplot(312)
plot(start:finish,(resrecord./max(resrecord(:,5:end)')')')
grid on
ylim([-5 1])
xlabel('Time [s]')
ylabel('Normalized weighted residuals')
legend('DQ','TBP')
subplot(313)
try
    plot(start:finish,diag([1 5e3])*(outcorrecord-[0; 273.15]),...
        start:finish,(outrecord-[0; 273.15])'*diag([1 5e3]),...
        [detectiontime detectiontime],[0 1.5e5],'--',...
        [detectiontime-500 detectiontime-500],[0 1.5e5],'--');
    legend('Reestimated DQ','Reestimated TBP*5e3',...
        'Faulty DQ','Faulty TBP*5e3','Estimated detection time','Parameter sampling')
    xlabel('Time [s]')
    ylabel('Outputs')
catch
    warning('No detection plot')
end
subplot(322)
fault = fault_sim.signals.values;
try
    if ~exist('fielddata')
        plot(start:finish,outrecord(2,:)-outcorrecord(2,:),...
            4000:finish,fault(4000:end-1)')
    else
        plot(start:finish,outrecord(2,:)-outcorrecord(2,:))
    end
catch
    warning('No fault plot')
end
xlabel('Time [s]')
ylabel('Fault estimation [K]')
legend('Estimated fault','True fault')
%% Residual
global var_resid FalseAlarmTime
if ~exist('fielddata')
    resid = resrecord(end,:);
else
    resid = [0 ew];
end
resid(1) = 0;
figure(11)
resid_normal = [0 resid(1,2:round(length(Y)/5))];
subplot(221)
autocorr(detrend(resid_normal))
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
% ------- Relative parameter standard deviations -------
outs = outrecord(:,2:round(length(Y)/5));
phis = phirecord(:,2:round(length(Y)/5));
if 1
%     figure(20)
%     colormap jet
    for it2 = 1:100:length(outs)-100
        paramvar1 = var(outs(1,it2:it2+100))./(phis(:,it2+100)*phis(:,it2+100)');
        paramvar2 = var(outs(2,it2:it2+100))./(phis(:,it2+100)*phis(:,it2+100)');
%         subplot(121), imagesc(paramvar1)
%         subplot(122), imagesc(paramvar2)
%         pause(0.1)
    end
end
diag(paramvar1)
diag(paramvar2)
var(paramrecord(:,length(outs)-200:length(outs)-100)')'
relparamstd1 = sqrt(diag(paramvar1))./paramrecord(1:length(phi),detectiontime-300-start)
relparamstd2 = sqrt(diag(paramvar2))./paramrecord(1+length(phi):2*length(phi),detectiontime-300-start)
% --------------- Detection ----------------------
figure(13)
subplot(311)
plot(start:finish,(resrecord./max(resrecord(:,5:end)')')')
ylim([-20 20])
grid on
ylabel('Residual')
xlabel('Time [s]')
subplot(334)
plot(start:finish,grecord(1,:)',start:finish,fdCUSUM.h*rectwin(length(Y)-start),'r--')
ylabel('g function')
xlabel('Time [s]')
title('CUSUM')
subplot(337)
plot(start:finish,grecord(1,:)'>fdCUSUM.h)
ylabel('Fault detection')
xlabel('Time [s]')
% GLR
subplot(335)
plot(start:finish,grecord(2,:)',start:finish,fdGLR.h*rectwin(length(Y)-start),'r--')
ylim([0 fdGLR.h*10])
ylabel('g function')
xlabel('Time [s]')
title('GLR')
subplot(338)
plot(start:finish,grecord(2,:)'>fdGLR.h,'LineWidth',2)
ylabel('Fault detection')
xlabel('Time [s]')
% EM
M = 500;
subplot(336)
plot(start:finish,grecord(3,:)',start:finish,fdEM.h*rectwin(length(Y)-start),'r--')
ylabel('g function')
xlabel('Time [s]')
title('Expectation Maximization (multidim.)')
subplot(339)
plot(start:finish,grecord(3,:)'>fdEM.h,'LineWidth',2)
ylabel('Fault detection')
xlabel('Time [s]')

if exist('fielddata')
    figure(15)
    plot(start+1:it,ew);
    title('Whitened residual')
end
