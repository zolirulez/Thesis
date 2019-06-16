clearvars
load resid_normal_chirp
% Input: [1; CRA; CRV; TA0; CRIT; THRd]
delay = 300;
Y = resid_normal;
U = Ur(2:length(resid_normal)+1);

% NOTE: for fielddata the investigated length is 1:7000, for 
% chirp it is 1:2000

y=Y(1:1000,1);
u=Ur(1:1000,4);
z=[y u];
Ze=dtrend(z);
y=Y(1001:2000,1);
u=Ur(1001:2000,4);
z=[y u];
Zt=dtrend(z);

nmax =3;
[nstruc1,losse1,losst1]=structbuilder(nmax,Ze,Zt);
[nstruc2,losse2,losst2]=structbuilder(nmax,Zt,Ze);

%% Exercise 1.1, 1.2
close all
N=length(Ze);
[AIC1 BIC1 FPE1 beststrt1] = findstruct(nstruc1,losse1,losst1,N,1);
paramanal(beststrt1,Zt,1);

% Backwards
[AIC2 BIC2 FPE2 beststrt2] = findstruct(nstruc2,losse2,losst2,N,2);
paramanal(beststrt2,Ze,2);