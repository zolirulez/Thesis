% This script generates the evaluation of the ARMAX filters on the
% residual, which for choosing an inverse whitening filter for fault
% detection

% clearvars
fielddata = 1;
if fielddata
    load resid_normal_fielddata2
else
    load resid_faultcontrol_ta0 
end
% Input: [1; CRA; CRV; TA0; CRIT; THRd]
Y = resrecord';
U = TA0; 

% NOTE: for fielddata the investigated length is 1:7000, for 
% chirp it is 1:2000
% for faultyestctrl it is 1:10000

% Inputs
if fielddata
    y=Y(1:2500,1);
    u=U(1:2500,1);
else
    y=Y(1:1500,1);
    u=U(1:1500,1);
end
z=[y u];
Ze=dtrend(z);
if fielddata
    y=Y(2501:5000,1);
    u=U(2501:5000,1);
else
    y=Y(1501:3000,1);
    u=U(1501:3000,1);
end
z=[y u];
Zt=dtrend(z);

% Building structures
nmax =3;
[nstruc1,losse1,losst1]=structbuilder(nmax,Ze,Zt);
[nstruc2,losse2,losst2]=structbuilder(nmax,Zt,Ze);


% -------- Find the best structure for every model order --------

% Forward
N=length(Ze);
[AIC1 BIC1 FPE1 beststrt1,h] = findstruct(nstruc1,losse1,losst1,N,1);
paramanal(beststrt1,Zt,1);

% Backwards
[AIC2 BIC2 FPE2 beststrt2,h] = findstruct(nstruc2,losse2,losst2,N,2);
paramanal(beststrt2,Ze,2);