clearvars
fielddata = 0;
if fielddata
%     load resid_normal_fielddata
    load resid_normal_fielddata2
else
%     load resid_normal_chirp
    load resid_faultcontrol_ta0 %resid_normal_faultyestctrl
end
% Input: [1; CRA; CRV; TA0; CRIT; THRd]
% delay = 75;
Y = resrecord(1,2000:5000)'; %Y = resid_normal;
U = TA0(2000:5000); %U = TA0(2:length(resid_normal)+1,1); %Ur(2:length(resid_normal)+1);

% NOTE: for fielddata the investigated length is 1:7000, for 
% chirp it is 1:2000
% for faultyestctrl it is 1:10000

if fielddata
    y=Y(1:3500,1);
    u=U(1:3500,1);
else
    y=Y(1:1500,1);
    u=U(1:1500,1);
end
z=[y u];
Ze=dtrend(z);
if fielddata
    y=Y(3501:7000,1);
    u=U(3501:7000,1);
else
    y=Y(1501:3000,1);
    u=U(1501:3000,1);
end
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