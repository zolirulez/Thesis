clearvars
load longsimulation_faulty
addpath('C:\Users\u375749\Documents\Thesis\Codes\ModelicaPlant\detection')
delta_h = record(8,:);
% Note: dont start it before 350, because of the delay operations
start = 2001;
resrecord = NaN(2,finish-start+1);
W = [5000; 6000; 1e3; 1e4; 1e3]*[1 10]; 
P = diag(max(W')')/10; 
for it = start:finish
    TA0 = U(it,10);
    THR = CoolProp.PropsSI('T','H',U(it,11),'P',Y(it,1),'CO2');
    CRIT = U(it-1,4)*VValues(3)*U(it-1,7);
    CRA = U(it-1,1);
    DV = CRA*DVValues(2);
    DmV = U(it,3)*KvValues*sqrt(U(it,6)*(Y(it,1) - Y(it,3)));
    DQ = DmV*(U(it,11) - hBP);
    DT = TBP - TA0;
    sigma = DQ/DT;
    phi = [1; U(it-1,1); TA0; CRIT; THR];
    res = [sigma delta_h(it-start+1)] - phi'*W;
    schur = 1 + phi'*P*phi;
    K = P*phi/schur ;
    W = W + K*res; 
    P = P - K*schur*K';
    resrecord(:,it-start+1) = res;
end

%% Residual
global var_resid FalseAlarmTime
mu0 = 0;
mu1 = 1e4;
resid = resrecord(2,:);
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
plot(start:finish,fault')
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
ylabel('g function')
xlabel('Time [s]')
title('GLR')
subplot(338)
plot(start:finish,fault')
ylabel('Fault detection')
xlabel('Time [s]')
% Moving square change detection
figure(7)
M = 500;
[g,fault,t,threshold] = MSCD(resid,M);
subplot(336)
plot(t+start,g',t+start,threshold*rectwin(length(resid)-M),'r--')
ylabel('Autocorrelation')
xlabel('Time [s]')
title('Moving Square Change Detector')
subplot(339)
plot(t+start,fault')
ylabel('Fault detection')
xlabel('Time [s]')