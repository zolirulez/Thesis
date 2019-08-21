% Script for plotting the estimation results of estimatortest_simulink

% For virtual state comparison
load TGC_faultcontrol.mat

% State estimation
t = start:Ts:finish;
tcw = start:Ts:it;
figure(1)
clf
subplot(221)
hold on
plot(t,record(1,:)/1e5)
plot(t,record(5,:)/1e5)
hold off
ylim([20 100])
xlabel('Time [s]')
ylabel('Pressure [bar]')
subplot(222)
hold on
plot(t,record(2,:)/1000)
plot(t,record(6,:)/1000)
plot(t,record(8,:)/1000)
hold off
ylim([-50 450])
xlabel('Time [s]')
ylabel('Enthalpy [kJ/kg]')
subplot(223)
hold on
plot(t,record(4,:)-273.15)
hold off
xlabel('Time [s]')
ylabel('Temperature [C]')
subplot(224)
hold on
plot(t,record(3,:))
plot(t,record(7,:))
hold off
ylim([100 600])
xlabel('Time [s]')
ylabel('Density [kg/m^3]')

% Linear filter performance
figure(2)
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
grid on
subplot(235)
plot(t,record(nx+nx+3,:)/1e3,t,record(nx+nx+6,:)/1e3)
legend('h_B_P','h_1')
xlabel('Time [s]')
ylabel('Innovation [kJ/kg]')
xlim([1000 length(Y)])
grid on
subplot(236)
plot(t,record(nx+nx+5,:)/1e3)
legend('h_R')
xlabel('Time [s]')
ylabel('Innovation [kJ/kg]')
xlim([1000 length(Y)])
grid on

% Earlier use: parameter estimation
% figure(3)
% plot(t,record(nx+nx+nx+ny+1:nx+nx+nx+ny+2,:))
% legend('s_0','k')
% xlabel('Time [s]')
% ylabel('Parameters')

figure(1)
% load hcw_chirp
subplot(221)
hold on
plot(t,recordf(1,:)/1e5,'--')
plot(t,recordf(5,:)/1e5,'--')
hold off
subplot(222)
hold on
plot(t,recordf(2,:)/1000,'--')
plot(t,recordf(6,:)/1000,'--')
plot(t,recordf(8,:)/1000,'--')
if ~exist('fielddata')
    plot(tcw,hcw(start:it)/1000,'--')
end
hold off
subplot(223)
hold on
plot(t,recordf(4,:)-273.15,'--')
if ~exist('fielddata')
    plot(tcw,Tcw(start:it)-273.15-3)
end
hold off
ylim([20 50])
xlabel('Time [s]')
ylabel('Temperature [C]')
subplot(224)
hold on
plot(t,recordf(3,:),'--')
plot(t,recordf(7,:),'--')
hold off

% Estimation versus true values plot
try
    figure(4)
    clf
    subplot(121)
    hold on
    plot(smooth(hcw(start+500:it)/1000,100),smooth(record(2,501:end)/1000,100))
    plot(smooth(hcw(start+500:it)/1000,100),smooth(recordf(2,501:end)/1000,100))
    plot(linspace(340,410,2),linspace(340,410,2))
    hold off
    xlabel('Theoretical approx. of enthalpy [kJ/kg]')
    ylabel('Estimation [kJ/kg]')
    legend('Original','Reestimated')
    title('Lissajous figure for enthalpy')
    xlim([340 390])
    ylim([340 390])
    subplot(122)
    hold on
    plot(smooth(Tcw(start+500:it)-273.15-3,100),smooth(record(4,501:end)-273.15-3,100))
    plot(smooth(Tcw(start+500:it)-273.15-3,100),smooth(recordf(4,501:end)-273.15-3,100))
    plot(linspace(25,50,2),linspace(25,50,2))
    hold off
    xlabel('Theoretical approx. of temperature [C]')
    ylabel('Estimation [C]')
    legend('Original','Reestimated')
    title('Lissajous figure for temperature')
    xlim([25 40])
    ylim([25 40])
    hold off
catch
    warning('No Lissajous figure')
end
    
% figure(3)
% hold on
% plot(t,recordf(nx+nx+nx+ny+1:nx+nx+nx+ny+2,:),'--')
% hold off
% legend('s_0','k')
% xlabel('Time [s]')
% ylabel('Parameters')

% Detection plot
figure(10)
try
    plot(0:finish,(outcorrecord(1:2,:)),...
        0:finish,(outrecord)',...
        [detectiontime detectiontime],[0 1.5e5],'--',...
        [detectiontime-500 detectiontime-500],[0 1.5e5],'--');
    legend({'Reestimated $\dot{Q}$','Reestimated $h_{BP}$',...
        'Faulty $\dot{Q}$','Faulty $h_{BP}$','Estimated detection time','Parameter sampling'},...
        'Interpreter','latex','Location','east')
    xlabel('Time [s]')
    ylabel('Outputs')
    xlim([1000 length(Y)])
catch
    warning('No detection plot')
end

% Residual
for it2 = 1:2
    figure(11)
    resid = resrecord(1,:);
    resid = resid(:);
    resid(1) = 0;
    resid_normal = [0 resid(2000:4000,1)'];
    subplot(221)
    autocorr(detrend(resid_normal))
    subplot(222)
    [pxx,f,pxxc] = periodogram(resid_normal,rectwin(length(resid_normal)),...
        length(resid_normal),1,'ConfidenceLevel',0.99);
    cpxx = cumsum(pxx)./sum(pxx);
    plot(f,cpxx,'b',f,f./max(f)+0.05,'r--',f+0.05*max(f),f./max(f),'r--')
    xlabel('Frequency')
    title('Cum. periodrogram with 99% conf. interval')
    subplot(223)
    probplot('normal',resid_normal)
    title('Prob. plot for normal dist.')
    subplot(224)
    var_resid1 = var(resid_normal);
    histfit(resid_normal,50)
    ylabel('Sample frequency')
    title('Histogram of data')
end

% % ------- Relative parameter standard deviations -------
outs = outrecord(:,2:round(length(Y)/5));
phis = phirecord(:,2:round(length(Y)/5));
for it2 = 1:100:length(outs)-100
    paramvar1 = var(outs(1,it2:it2+100))./(phis(:,it2+100)*phis(:,it2+100)');
    paramvar2 = var(outs(2,it2:it2+100))./(phis(:,it2+100)*phis(:,it2+100)');
end
diag(paramvar1)
diag(paramvar2)
var(paramrecord(:,length(outs)-200:length(outs)-100)')'
relparamstd1 = sqrt(diag(paramvar1))./paramrecord(1:length(phi),detectiontime-300-start)
relparamstd2 = sqrt(diag(paramvar2))./paramrecord(1+length(phi):2*length(phi),detectiontime-300-start)

% --------------- Detection ----------------------
figure(13)
clf
subplot(321)
hold on
plot(start+1:it,resrecord(1,start+1:it))
hold off
ylabel('Whitened residual')
xlabel('Time [s]')
xlim([2000 length(Y)])
ylabel('Residual of $\dot{Q}$','Interpreter','latex')
xlabel('Time [s]')
legend('Original','Whitened')
subplot(322)
plot(start+1:it,resrecord(2,start+1:it))
ylabel('Whitened residual')
xlabel('Time [s]')
xlim([2000 length(Y)])
ylabel('Residual of h_B_P')
xlabel('Time [s]')
subplot(334)
plot(2000:finish,grecord(1,2001:end)',0:finish,fdCUSUM.h*rectwin(length(Y)),'r--')
ylabel('g function')
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
ylabel('g function')
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
ylabel('g function')
xlabel('Time [s]')
title('Expectation Maximization (multidim.)')
xlim([2000 length(Y)])
subplot(339)
plot(2000:finish,grecord(3,2001:end)'>fdEM.h,'LineWidth',2)
ylabel('Fault detection')
xlabel('Time [s]')
xlim([2000 length(Y)])

% Measurement reconstruction
figure(5)
clf
hold on
plot(1000:finish,Y(1001:end,2)/1000)
plot(1000:finish,outcorrecord(4,1001:end)/1000,'LineWidth',2)
plot(1000:finish,outcorrecord(2,1001:end)/1000)
plot(1000:finish,outcorrecord(3,1001:end)/1000)
plot([detectiontime detectiontime],[0 1.5e5],'--',...
    [detectiontime-500 detectiontime-500],[0 1.5e5],'--',...
    [5000 5000],[0 1.5e5],'--')
hold off
ylabel('h_B_P [kJ/kg]')
xlabel('Time [s]')
handle = legend({'$h_{BPm}$','$h_{BP}$','$\hat{h}_{BP}$',...
    '$h_{HR}-\hat{Q}/\dot{m}_{HR}$','Detection time','Parameter sampling',...
    'Fault occurence'},'Location','southeast'...
    ,'Interpreter','Latex');
ylim([230 330])
xlim([1000 length(Y)])


