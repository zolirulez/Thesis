clearvars
addpath('C:\Users\u375749\Documents\Thesis\Codes\ModelicaPlant\bestsimulations')
addpath('C:\Users\u375749\Documents\Thesis\Codes\ModelicaPlant\objects_inits')
addpath('C:\Users\u375749\Documents\Thesis\Codes\Measurements')
TsParam = 1000;

dprime = NaN(3,5);
falsealarm = NaN(3,3);
misseddetection = NaN(3,3);

h = figure(1);
set(h, 'Position',  [100, 100, 100+1400, 100+300])
clf
subplot(131)
simulationstring = 'faultignore_ta0';
estimatortest_simulink
hold on
[FPR,TPR,T,AUCRES1] = perfcurve([zeros(5001,1); ones(9000,1); zeros(6000,1)],...
    resrecord(1,:)',1);
plot(FPR,TPR)
[FPR,TPR,T,AUCRES2] = perfcurve([zeros(5001,1); ones(9000,1); zeros(6000,1)],...
    -resrecord(2,:)',1);
plot(FPR,TPR)
[FPR,TPR,T,AUCCUSUM] = perfcurve([zeros(5001,1); ones(9000,1); zeros(6000,1)],...
    grecord(1,:)',1);
plot(FPR,TPR)
falsealarmCUSUM = FPR(find(T<fdCUSUM.h,1,'first'));
misseddetectionCUSUM = 1-TPR(find(T<fdCUSUM.h,1,'first'));
[FPR,TPR,T,AUCGLR] = perfcurve([zeros(5001,1); ones(9000,1); zeros(6000,1)],...
    grecord(2,:)',1);
plot(FPR,TPR)
falsealarmGLR = FPR(find(T<fdGLR.h,1,'first'));
misseddetectionGLR = 1-TPR(find(T<fdGLR.h,1,'first'));
[FPR,TPR,T,AUCEM] = perfcurve([zeros(5001,1); ones(9000,1); zeros(6000,1)],...
    grecord(3,:)',1);
plot(FPR,TPR)
falsealarmEM = FPR(find(T<0.5,1,'first'));
misseddetectionEM = 1-TPR(find(T<0.5,1,'first'));
plot([0 1],[0 1],'--')
hold off
legend({['$\dot{\tilde{Q}}$, AUC:' num2str(AUCRES1)],...
    ['$\tilde{h}_{BP}$, AUC:' num2str(AUCRES2)],...
    ['CUSUM, AUC: ' num2str(AUCCUSUM)],['GLR, AUC: ' num2str(AUCGLR)],...
    ['EM, AUC: ' num2str(AUCEM)],'Random classifier, AUC: 0.5'},...
    'Location','southeast','Interpreter','latex')
xlabel('False alarm rate')
ylabel('Hit rate')
title('Ignored fault, simulation')
dprime(1,:) = sqrt(2)*[norminv(AUCRES1) norminv(AUCRES2)...
    norminv(AUCCUSUM) norminv(AUCGLR) norminv(AUCEM)];
falsealarm(1,:) = [falsealarmCUSUM falsealarmGLR falsealarmEM];
misseddetection(1,:) = [misseddetectionCUSUM misseddetectionGLR misseddetectionEM];

subplot(132)
simulationstring = 'faultcontrol_ta0';
estimatortest_simulink
hold on
[FPR,TPR,T,AUCRES1] = perfcurve([zeros(5001,1); ones(9000,1); zeros(6000,1)],...
    resrecord(1,:)',1);
plot(FPR,TPR)
[FPR,TPR,T,AUCRES2] = perfcurve([zeros(5001,1); ones(9000,1); zeros(6000,1)],...
    -resrecord(2,:)',1);
plot(FPR,TPR)
[FPR,TPR,T,AUCCUSUM] = perfcurve([zeros(5001,1); ones(9000,1); zeros(6000,1)],...
    grecord(1,:)',1);
plot(FPR,TPR)
falsealarmCUSUM = FPR(find(T<fdCUSUM.h,1,'first'));
misseddetectionCUSUM = 1-TPR(find(T<fdCUSUM.h,1,'first'));
[FPR,TPR,T,AUCGLR] = perfcurve([zeros(5001,1); ones(9000,1); zeros(6000,1)],...
    grecord(2,:)',1);
plot(FPR,TPR)
falsealarmGLR = FPR(find(T<fdGLR.h,1,'first'));
misseddetectionGLR = 1-TPR(find(T<fdGLR.h,1,'first'));
[FPR,TPR,T,AUCEM] = perfcurve([zeros(5001,1); ones(9000,1); zeros(6000,1)],...
    grecord(3,:)',1);
plot(FPR,TPR)
falsealarmEM = FPR(find(T<0.5,1,'first'));
misseddetectionEM = 1-TPR(find(T<0.5,1,'first'));
plot([0 1],[0 1],'--')
hold off
legend({['$\dot{\tilde{Q}}$, AUC:' num2str(AUCRES1)],...
    ['$\tilde{h}_{BP}$, AUC:' num2str(AUCRES2)],...
    ['CUSUM, AUC: ' num2str(AUCCUSUM)],['GLR, AUC: ' num2str(AUCGLR)],...
    ['EM, AUC: ' num2str(AUCEM)],'Random classifier, AUC: 0.5'},...
    'Location','southeast','Interpreter','latex')
xlabel('False alarm rate')
ylabel('Hit rate')
title('Fault handling, simulation')
dprime(2,:) = sqrt(2)*[norminv(AUCRES1) norminv(AUCRES2)...
    norminv(AUCCUSUM) norminv(AUCGLR) norminv(AUCEM)];
falsealarm(2,:) = [falsealarmCUSUM falsealarmGLR falsealarmEM];
misseddetection(2,:) = [misseddetectionCUSUM misseddetectionGLR misseddetectionEM];

subplot(133)
main_meas
hold on
[FPR,TPR,T,AUCRES1] = perfcurve([zeros(7747-1001,1); ones(12367-7747,1); zeros(length(Y)-12367,1)],...
    resrecord(1,:)',1);
plot(FPR,TPR)
[FPR,TPR,T,AUCRES2] = perfcurve([zeros(7747-1001,1); ones(12367-7747,1); zeros(length(Y)-12367,1)],...
    -resrecord(2,:)',1);
plot(FPR,TPR)
[FPR,TPR,T,AUCCUSUM] = perfcurve([zeros(7747-1001,1); ones(12367-7747,1); zeros(length(Y)-12367,1)],...
    grecord(1,:)',1);
plot(FPR,TPR)
falsealarmCUSUM = FPR(find(T<fdCUSUM.h,1,'first'));
misseddetectionCUSUM = 1-TPR(find(T<fdCUSUM.h,1,'first'));
[FPR,TPR,T,AUCGLR] = perfcurve([zeros(7747-1001,1); ones(12367-7747,1); zeros(length(Y)-12367,1)],...
    grecord(2,:)',1);
plot(FPR,TPR)
falsealarmGLR = FPR(find(T<fdGLR.h,1,'first'));
misseddetectionGLR = 1-TPR(find(T<fdGLR.h,1,'first'));
[FPR,TPR,T,AUCEM] = perfcurve([zeros(7747-1001,1); ones(12367-7747,1); zeros(length(Y)-12367,1)],...
    grecord(3,:)',1);
plot(FPR,TPR)
falsealarmEM = FPR(find(T<0.5,1,'first'));
misseddetectionEM = 1-TPR(find(T<0.5,1,'first'));
plot([0 1],[0 1],'--')
hold off
legend({['$\dot{\tilde{Q}}$, AUC:' num2str(AUCRES1)],...
    ['$\tilde{h}_{BP}$, AUC:' num2str(AUCRES2)],...
    ['CUSUM, AUC: ' num2str(AUCCUSUM)],['GLR, AUC: ' num2str(AUCGLR)],...
    ['EM, AUC: ' num2str(AUCEM)],'Random classifier, AUC: 0.5'},...
    'Location','southeast','Interpreter','latex')
xlabel('False alarm rate')
ylabel('Hit rate')
title('Field data (ignored fault)')
dprime(3,:) = sqrt(2)*[norminv(AUCRES1) norminv(AUCRES2)...
    norminv(AUCCUSUM) norminv(AUCGLR) norminv(AUCEM)];
falsealarm(3,:) = [falsealarmCUSUM falsealarmGLR falsealarmEM];
misseddetection(3,:) = [misseddetectionCUSUM misseddetectionGLR misseddetectionEM];

saveas(h,'roc.png')
save('dprime','dprime')
dprime