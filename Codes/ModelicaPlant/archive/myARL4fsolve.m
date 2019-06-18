function [resid] = myARL4fsolve(mu1_treshold_vector)

% Input demux
mu1 = mu1_treshold_vector(1);
threshold = mu1_treshold_vector(2);

% Given parameters, not to be changed
mu0 = 0;
global var_resid FalseAlarmTime

% myARL function
h = threshold;
sigs = sqrt((mu1-mu0)^2/var_resid);
mus = sqrt((mu1-mu0)^2/2/var_resid);
L = @(mus,sigs,h)(sigs^2/2/mus^2*(exp(-2*(mus*h/sigs^2+1.166*mus/sigs))-1+2*(mus*h/sigs^2+1.166*mus/sigs)));
Tdetect = L(mus^2,sigs^2,h);
Tfalse = L(-mus^2,sigs^2,h);

% Checking how close we are to the solution
detectresid = Tdetect;
falseresid = Tfalse - FalseAlarmTime;
resid = [detectresid; falseresid];