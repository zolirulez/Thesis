clearvars
addpath('C:\Users\u375749\Documents\Thesis\Codes\ModelicaPlant\bestsimulations')
addpath('C:\Users\u375749\Documents\Thesis\Codes\ModelicaPlant\objects_inits')
addpath('C:\Users\u375749\Documents\Thesis\Codes\Measurements')
TsParam = 1000;

dprime = NaN(3,3);

h = figure(1);
set(h, 'Position',  [100, 100, 100+1200, 100+300])
clf
subplot(131)
simulationstring = 'faultignore_ta0';
estimatortest_simulink
hold on
[FPR,TPR,T,AUCCUSUM] = perfcurve([zeros(5001,1); ones(9000,1); zeros(6000,1)],...
    grecord(1,:)',1);
plot(FPR,TPR)
[FPR,TPR,T,AUCGLR] = perfcurve([zeros(5001,1); ones(9000,1); zeros(6000,1)],...
    grecord(2,:)',1);
plot(FPR,TPR)
[FPR,TPR,T,AUCEM] = perfcurve([zeros(5001,1); ones(9000,1); zeros(6000,1)],...
    grecord(3,:)',1);
plot(FPR,TPR)
plot([0 1],[0 1],'--')
hold off
legend({['CUSUM, AUC: ' num2str(AUCCUSUM)],['GLR, AUC: ' num2str(AUCGLR)],...
    ['EM, AUC: ' num2str(AUCEM)],'Random classifier ROC'},'Location','southeast')
xlabel('False positive rate')
ylabel('True positive rate')
title('Ignored fault, simulation')
dprime(1,:) = sqrt(2)*[norminv(AUCCUSUM) norminv(AUCGLR) norminv(AUCEM)];

subplot(132)
simulationstring = 'faultcontrol_ta0';
estimatortest_simulink
hold on
[FPR,TPR,T,AUCCUSUM] = perfcurve([zeros(5001,1); ones(9000,1); zeros(6000,1)],...
    grecord(1,:)',1);
plot(FPR,TPR)
[FPR,TPR,T,AUCGLR] = perfcurve([zeros(5001,1); ones(9000,1); zeros(6000,1)],...
    grecord(2,:)',1);
plot(FPR,TPR)
[FPR,TPR,T,AUCEM] = perfcurve([zeros(5001,1); ones(9000,1); zeros(6000,1)],...
    grecord(3,:)',1);
plot(FPR,TPR)
plot([0 1],[0 1],'--')
hold off
legend({['CUSUM, AUC: ' num2str(AUCCUSUM)],['GLR, AUC: ' num2str(AUCGLR)],...
    ['EM, AUC: ' num2str(AUCEM)],'Random classifier ROC'},'Location','southeast')
xlabel('False positive rate')
ylabel('True positive rate')
title('Fault operation, simulation')
dprime(2,:) = sqrt(2)*[norminv(AUCCUSUM) norminv(AUCGLR) norminv(AUCEM)];

subplot(133)
main_meas
hold on
[FPR,TPR,T,AUCCUSUM] = perfcurve([zeros(7747-1001,1); ones(12367-7747,1); zeros(length(Y)-12367,1)],...
    grecord(1,:)',1);
plot(FPR,TPR)
[FPR,TPR,T,AUCGLR] = perfcurve([zeros(7747-1001,1); ones(12367-7747,1); zeros(length(Y)-12367,1)],...
    grecord(2,:)',1);
plot(FPR,TPR)
[FPR,TPR,T,AUCEM] = perfcurve([zeros(7747-1001,1); ones(12367-7747,1); zeros(length(Y)-12367,1)],...
    grecord(3,:)',1);
plot(FPR,TPR)
plot([0 1],[0 1],'--')
hold off
legend({['CUSUM, AUC: ' num2str(AUCCUSUM)],['GLR, AUC: ' num2str(AUCGLR)],...
    ['EM, AUC: ' num2str(AUCEM)],'Random classifier ROC'},'Location','southeast')
xlabel('False positive rate')
ylabel('True positive rate')
title('Field data (ignored fault)')
dprime(3,:) = sqrt(2)*[norminv(AUCCUSUM) norminv(AUCGLR) norminv(AUCEM)];

saveas(h,'roc.png')
save('dprime','dprime')
dprime