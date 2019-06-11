clearvars
load uy_sim_faulty%longsimulation_faulty
load fault_sim
addpath('C:\Users\u375749\Documents\Thesis\Codes\ModelicaPlant\detection')
nu = 12;
ny = 5;
U = uy_sim.signals.values(:,1:nu);
Y = uy_sim.signals.values(:,nu+1:nu+ny);
% Note: dont start it before 350, because of the delay operations
start = 2001;
finish = 10000;
rlsInitial.t = [5000; 6000; 1e3; 1e4; 1e3]*[1 1e-3]; % 10
rlsInitial.P = diag(max(rlsInitial.t')')/10; 
rls = RecursiveLeastSquares;
lambda = 1-1e-5;
weight = eye(2);
rls.initialize(lambda,weight,rlsInitial);

% Fault Detector
fdGLR = FaultDetector;
mean.m0 = 0;
variance = 0.25; %1.8614e+06;
param.WindowLength = 500;
param.InitialGuess = [-1,1];
% param.Threshold = 500;
param.FalseAlarmProbability = 0*1e-20;
method = 'GLR';
fdGLR.initialize(mean,variance,method,param);
clear param
fdCUSUM = FaultDetector;
mean.m1 = -1;%1.0270e+04;
%param.Threshold = 51.1729;
param.FalseAlarmTime = 1e20;
param.InitialGuess = [mean.m1,1];
method = 'CUSUM';
fdCUSUM.initialize(mean,variance,method,param);
clear param
fdEM = FaultDetector;
method = 'EM';
param.MeanDensityRatio = 0.1;
param.SamplingTime = 1;
param.ResponsibilityTimeConstant = 100;
mean.m0 = zeros(2,1);
variance = diag([5e6 5]);%diag([4.25e6 0.25]);
fdEM.initialize(mean,variance,method,param);

% Measurement correction
% It seems that the expectable detection is around 60 s.
% Choosing safety parameter 5, we take W 300 s before.
detectiontime = 4060;
savetime = detectiontime - 300;
faultOperation = 0;

% Recording
resrecord = NaN(2,finish-start+1);
paramrecord = NaN(size(rls.t,1)*size(rls.t,2),finish-start+1);
outrecord = NaN(2,finish-start+1);
outcorrecord = NaN(2,finish-start+1);
% outcorrecord = [];
grecord = NaN(3,finish-start+1);
faultrecord = NaN(3,finish-start+1);
for it = start:finish
    TA0 = U(it,10);
    THR = CoolProp.PropsSI('T','H',U(it-300,11),'P',Y(it-300,1),'CO2');
    hBP = Y(it,2);
    TBP = CoolProp.PropsSI('T','H',hBP,'P',Y(it,1),'CO2');
    % DmG = U(it-1,4)*VValues(3)*U(it-1,7);
    CRIT = U(it-300,4);
    CRA = U(it-300,1);
%    DV = CRA*DVValues(2);
    DmV = U(it,3)*2.4e-5*sqrt(U(it,6)*(Y(it,1) - Y(it,3)));
    DQ = DmV*(U(it-300,11) - hBP);
    DT = TBP - TA0;
    sigma = DQ/DT;
    phi = [1; CRA; TA0; CRIT; THR];
    out = [DQ DT];
    rls.regression(phi,out);
    if it == start
        rls.e = [0 0];
    end
    [~,fault1] = fdCUSUM.CUSUM(rls.e(2));
    [~,fault2] = fdGLR.GLR(rls.e(2));
    [~,fault3] = fdEM.EM(rls.e');
    % Measurement correction: based on EM fault detector
    if it > start + 20
        if sum(faultrecord(3,it-19-start:it-10-start)) < 2 && all(faultrecord(3,it-9-start:it-start))
            Wsave = reshape(paramrecord(:,it-start),length(phi),length(out));
            detectiontime = it;
            faultOperation = 1;
        end
        if sum(faultrecord(3,it-19-start:it-10-start)) > 8 && ~any(faultrecord(3,it-9-start:it-start))
            faultOperation = 0;
        end
    end
    if faultOperation
%         Wp = W'/(W*W'); % right pseudoinverse
%         outcor = out*Wp*Wsave;
        outcor = phi'*Wsave;
        outcorrecord(:,it-start+1) = outcor';
    end
    resrecord(:,it-start+1) = rls.e;
    paramrecord(:,it-start+1) = reshape(rls.t,size(rls.t,1)*size(rls.t,2),1);
    outrecord(:,it-start+1) = out';
    grecord(:,it-start+1) = [fdCUSUM.g; fdGLR.g; fdEM.g];
    faultrecord(:,it-start+1) = [fault1; fault2; fault3];
end
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
legend('DQ','DT - used now')
subplot(313)
try
plot(detectiontime:length(outcorrecord)+detectiontime-1,diag([1 5e3])*outcorrecord',...
    start:finish,outrecord'*diag([1 5e3]),...
    [detectiontime detectiontime],[0 1.5e5],'--',...
    [detectiontime-300 detectiontime-300],[0 1.5e5],'--');
catch
    warning('No detection plot')
end
legend('Reestimated DQ','Reestimated DT*5e3',...
    'Faulty DQ','Faulty DT*5e3','Estimated detection time','Parameter sampling')
xlabel('Time [s]')
ylabel('Outputs')
subplot(322)
fault = fault_sim.signals.values;
try
    plot(start:finish,outrecord(2,:)-outcorrecord(2,:),...
        4000:finish,fault(4000:end-1)')
catch
    warning('No fault plot')
end
xlabel('Time [s]')
ylabel('Fault estimation [K]')
legend('Estimated fault','True fault')
%% Residual
global var_resid
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
% 3 month = 3*30 days = 3*30*24 hours = 3*30*24*3600 seconds
% mu0 = 0;
% mu1 = 1e4;
% FalseAlarmTime = 10;
% initialguess = [mu1; 1];
% mu1_threshold_vector = fsolve('myARL4fsolve',initialguess);
% mu1 = mu1_threshold_vector(1)
% threshold = mu1_threshold_vector(2)
% [Tdetect, Tfalse] = myARL(mu1,mu0,var_resid,threshold);
% [g,fault,threshold] = myCusum(mu1,mu0,var_resid,resid,threshold);
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
% threshold = 500;
% M = 500;
% [Pmissed, Pfalse] = GLR_design(M, mu1,mu0,sqrt(var_resid),threshold,0)
% [g,fault] = myGLR(mu0,var_resid,resid,M,threshold);
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
% [g,fault,t,threshold] = MSCD(resid,M);
subplot(336)
plot(start:finish,grecord(3,:)',start:finish,fdEM.h*rectwin(length(resid)),'r--')
ylabel('g function')
xlabel('Time [s]')
title('Expectation Maximization (multidim.)')
subplot(339)
plot(start:finish,grecord(3,:)'>fdEM.h,'LineWidth',2)
ylabel('Fault detection')
xlabel('Time [s]')

