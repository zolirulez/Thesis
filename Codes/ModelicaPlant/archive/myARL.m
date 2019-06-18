function [Tdetect,Tfalse] = myARL(mu1,mu0,var,threshold)

h = threshold;
sigs = sqrt((mu1-mu0)^2/var); 
mus = sqrt((mu1-mu0)^2/2/var);
L = @(mus,sigs,h)(sigs^2/2/mus^2*(exp(-2*(mus*h/sigs^2+1.166*mus/sigs))-1+2*(mus*h/sigs^2+1.166*mus/sigs)));
Tdetect = L(mus,sigs,h);
Tfalse = L(-mus,sigs,h);