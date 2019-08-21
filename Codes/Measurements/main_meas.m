% For evaluation, the handle and the statistics
clearvars -except h dprime fdCUSUM fdGLR falsealarm misseddetection
% Reading data
load fielddata
addpath('C:\Users\u375749\Documents\Thesis\Codes\ModelicaPlant')
% Estimation and fault detection
estimatortest
disp(['Detection time ' num2str(detectiontime-7748)])
disp(['Fault switch off detection time ' num2str(switchofftime-12367)])