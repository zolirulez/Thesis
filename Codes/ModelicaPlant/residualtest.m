clearvars
load longsimulation_faulty
load fault_sim
addpath('C:\Users\u375749\Documents\Thesis\Codes\ModelicaPlant\detection')
delta_h = record(8,:);
% Note: dont start it before 350, because of the delay operations
start = 2001;
W = [5000; 6000; 1e3; 1e4; 1e3]*[1 1e-3 10]; 
P = diag(max(W')')/10; 

% Measurement correction
% It seems that the expectable detection is around 60 s.
% Choosing safety parameter 5, we take W 300 s before.
detectiontime = 4060;
savetime = detectiontime - 300;

% Recording
resrecord = NaN(3,finish-start+1);
paramrecord = NaN(size(W,1)*size(W,2),finish-start+1);
outrecord = NaN(3,finish-start+1);
outcorrecord = NaN(3,finish-detectiontime+1);
%Paramception
% load W2
% W2save = W2;
% W2 = reshape(W2,5,10);
% P2 = 1e10*eye(5);
% param2record = NaN(size(W2,1)*size(W2,2),finish-start+1);
for it = start:finish
    TA0 = U(it,10);
    THR = CoolProp.PropsSI('T','H',U(it,11),'P',Y(it,1),'CO2');
    hBP = Y(it,2);
    TBP = CoolProp.PropsSI('T','H',hBP,'P',Y(it,1),'CO2');
    CRIT = U(it-1,4)*VValues(3)*U(it-1,7);
    CRA = U(it-1,1);
    DV = CRA*DVValues(2);
    DmV = U(it,3)*KvValues*sqrt(U(it,6)*(Y(it,1) - Y(it,3)));
    DQ = DmV*(U(it,11) - hBP);
    DT = max(1,TBP - TA0);
    sigma = DQ/DT;
    phi = [1; CRA; TA0; CRIT; THR];
    out = [DQ DT delta_h(it-start+1)];
    res = out - phi'*W;
    res = res*diag([1 1 1]);
    schur = 1 + phi'*P*phi;
    K = P*phi/schur ;
    W = W + K*res; 
    P = P - K*schur*K';
%     % Paramception
%     out2 = reshape(W,1,size(W,1)*size(W,2));
%     res2 = out2 - phi'*W2;
%     schur2 = 1 + phi'*P2*phi;
%     K2 = P2*phi/schur2;
%     W2 = W2 + K2*res2;
%     P2 = P2 - K2*schur2*K2';
%     param2record(:,it-start+1) = reshape(W2,size(W2,1)*size(W2,2),1);
    if it == 3750
        disp('')
    end
    % Measurement correction
    if it == detectiontime
        Wsave = reshape(paramrecord(:,savetime-start),length(phi),length(out));
%         Wsave = phi'*W2save;
    end
    if it >= detectiontime
%         Wp = W'/(W*W'); % right pseudoinverse
%         outcor = out*Wp*Wsave;
        outcor = phi'*Wsave; %reshape(Wsave,size(W,1),size(W,2));
        outcorrecord(:,it-detectiontime+1) = outcor;
    end
    resrecord(:,it-start+1) = res;
    paramrecord(:,it-start+1) = reshape(W,size(W,1)*size(W,2),1);
    outrecord(:,it-start+1) = out';
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
legend('DQ','DT','delta_h - used now')
subplot(313)
plot(start+2001:length(outcorrecord)+start+2000,outcorrecord'*diag([1 5e3 1]),...
    start:finish,outrecord'*diag([1 5e3 1]),...
    [detectiontime detectiontime],[0 1.5e5],'--',...
    [detectiontime-300 detectiontime-300],[0 1.5e5],'--');
legend('Reestimated DQ','Reestimated DT*5e3','Reestimated delta_h',...
    'Faulty DQ','Faulty DT*5e3','Faulty delta_h','Estimated detection time','Parameter sampling')
xlabel('Time [s]')
ylabel('Outputs')
subplot(322)
fault = fault_sim.signals.values;
plot(start+2001:length(outcorrecord)+start+2000,outrecord(2,detectiontime-start+1:end)-outcorrecord(2,:),...
    start+2001:length(outcorrecord)+start+2000,fault(detectiontime+1:end)')
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
% 3 month = 3*30 days = 3*30*24 hours = 3*30*24*3600 seconds
mu0 = 0;
mu1 = 1e4;
FalseAlarmTime = 10;
initialguess = [mu1; 1];
mu1_threshold_vector = fsolve('myARL4fsolve',initialguess);
mu1 = mu1_threshold_vector(1)
threshold = mu1_threshold_vector(2)
[Tdetect, Tfalse] = myARL(mu1,mu0,var_resid,threshold)
[g,fault,threshold] = myCusum(mu1,mu0,var_resid,resid,threshold);
figure(7)
subplot(311)
plot(start:finish,resid')
grid on
ylabel('Residual')
xlabel('Time [s]')
ylim([-2e4 2e4])
subplot(334)
plot(start:finish,g',start:finish,threshold*rectwin(length(resid)),'r--')
ylabel('g function')
xlabel('Time [s]')
title('CUSUM')
subplot(337)
plot(start:finish,fault','LineWidth',2)
ylabel('Fault detection')
xlabel('Time [s]')
threshold = 500;
M = 500;
[Pmissed, Pfalse] = GLR_design(M, mu1,mu0,sqrt(var_resid),threshold,0)
[g,fault] = myGLR(mu0,var_resid,resid,M,threshold);
% GLR
figure(7)
subplot(335)
plot(start:finish,g',start:finish,threshold*rectwin(length(resid)),'r--')
ylim([0 threshold*10])
ylabel('g function')
xlabel('Time [s]')
title('GLR')
subplot(338)
plot(start:finish,fault','LineWidth',2)
ylabel('Fault detection')
xlabel('Time [s]')
% Moving square change detection
figure(7)
M = 500;
[g,fault,t,threshold] = MSCD(resid,M);
subplot(336)
plot(t+start,g',t+start,threshold*rectwin(length(resid)-M),'r--')
ylim([0 threshold*10])
ylabel('Autocorrelation')
xlabel('Time [s]')
title('Moving Square Change Detector')
subplot(339)
plot(t+start,fault','LineWidth',2)
ylabel('Fault detection')
xlabel('Time [s]')

