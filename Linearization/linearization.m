clearvars
close all
syms KvHP KvG                           % Kv value of valve
syms TauHP TauG TauVA TauBP TauBC TauR TauTA TauIT  % Time constants
syms R                                  % Hydraulic resistance
syms Vc VR VIT                          % Volumes of cells, reseiver and piston
syms DdDp1 DdDh1 DdDp2 DdDh2 DdDpR DdDhR DdDpBP DdDhBP
syms eS                                 % Isentropic efficiency
syms p1 p2 pR pMT                       % Pressures
syms h1 h2 hL hG hR hMT hHR hC hF       % Enthalpies
syms d1 d2 dR dG dA dBP                 % Densities
syms CRHP CRG BPR BP                    % Capacity or division ratios
syms fIT fA                             % Frequencies
syms DmHP DmIT DmG DmMT Dm21            % Mass flow rates
syms DVA                                % Volume flow rate of air
syms s s0 k                             % Convection parameters
syms cp                                 % Specific heat capacity of air
syms TA1 TA2 T1 T2 TBP TA0              % Temperatures
syms DTDp1 DTDh1 DTDp2 DTDh2 DTDpBP DTDhBP
syms DhDpIT                             % Isentropic curve delta ratio
syms DQC DQF                            % Heat flow rates
syms MxDVA                              % Maximal volume flow of fan
syms delta_hHR delta_h2
% ----------------------- STATIC EQUATIONS --------------------------------
% Note the parameters for identification TODO hL dBP dIT
% Boundary condition
DmC = DQC/(hC-hL);
DmF = DQF/(hF-hL);
DmL = DmC + DmF;
% Compressor
hIT = hG + DhDpIT*(p1 - pR)/eS;
% Heat recovery
DmHR = DmMT + DmIT;
hHR = (hMT*DmMT+hIT*DmIT)/DmHR - delta_hHR;
% Gas cooler, heat transfer
s = s0 + k*DVA;
w = dA*DVA*cp/s;
DQ1 = (TA1 - T1)*s;
DQ2 = (TA2 - T2)*s;
% Gas cooler, fluid side
Joining1 = [-1 d1; DdDp1 DdDh1];
Joining2 = [-1 d2; DdDp2 DdDh2];
DmBP = DmHP*BP;
Dm1 = DmHR - DmBP;
Dm2 = DmHP*(1-BP);
DPsi1 = Dm1*hHR-Dm21*h1 + DQ1;
DPsi2 = Dm21*h1-Dm2*h2 + DQ2;
% ByPass valve
hBP = (Dm2*h2 + DmBP*hHR)/DmHP;
% High pressure valve
% Receiver
JointR = [-1 dR; DdDpR DdDhR];
DPsiR = 1/VR*(DmHP*hBP - DmL*hL - DmG*hG);
% Receiver Valve
% Fan
% ---------------------- DYNAMIC EQUATIONS --------------------------------
% Boundary condition
DDmMT = 1/TauBC*(-DmMT + DmL);
% Heat Recovery (no differentials)
% Gas cooler
DDm21 = 1/TauR*(-Dm21 + 1/R*sqrt(d1*(p1-p2)));
Dd1 = 1/Vc*(Dm1-Dm21);
Dd2 = 1/Vc*(Dm21-DmHP*(1-BP));
Dph1 = Joining1\[DPsi1/Vc; Dd1];
Dph2 = Joining2\[DPsi2/Vc; Dd2];
DT1 = [DTDp1 DTDh1]*Dph1;
DT2 = [DTDp2 DTDh2]*Dph2;
DTA1 = 1/TauTA*(-TA1 + 1/(w+1)*T1 + 1/(1/w+1)*TA0);
DTA2 = 1/TauTA*(-TA2 + 1/(w+1)*T2 + 1/(1/w+1)*TA1);
% ByPass Valve
DBP = 1/TauBP*(-BP + BPR);
% High pressure valve
DDmHP = 1/TauHP*(-DmHP + CRHP*KvHP*sqrt(dBP*(p2 - pR)));
% Receiver
DdR = 1/VR*(DmHP-DmL-DmG);
DphR = JointR\[DPsiR; DdR];
% Receiver valve
DDmG = 1/TauG*(-DmG + CRG*KvG*sqrt(dG*(pR - pMT)));
% Compressor
DDmIT = 1/TauIT*(-DmIT + dG*VIT*fIT);
% Fan (A refers to air)
DDVA = 1/TauVA*(-DVA + dA*MxDVA*fA);
% ------------------------ LINEARIZATION ----------------------------------
% Augmentation
Dx = [DDVA; Dph1; Dd1; DT1; Dph2; Dd2; DT2; DDmHP; DphR; DdR; DDmG; DphC; DdC; DdMT; DDmMT];
x = [DVA; p1; h1; d1; T1; p2; h2; d2; T2; DmHP; pR; hR; dR; DmG; pC; hC; dC; dMT; DmMT];
u = [fA; CRHP; CRG; fMT];
d = [dA; TA0; DmMT; hL; hG; dG; hCi; hF; DQC; DQF];
% Jacobians
A = jacobian(Dx,x);
B = jacobian(Dx,u);
G = jacobian(Dx,d);
% Substitutions check the values TODO
A = subs(A,{p1 p2 pR pC},num2cell([86 84.8 38 30]*10^5));
A = subs(A,{h1 h2 hMT hL hG hR hC hCi hF},num2cell([450 350 525 212 430 300 460 460 460]*10^3));
A = subs(A,{d1 d2 dR dG dC dMT dA},num2cell([280 700 220 100 80 80 1.25]));
A = subs(A,{TA0 TA1 TA2 T1 T2},num2cell([30 35 40 90 50]+273.15));
A = subs(A,{DmHP DmMT DmG DmL DmCi DmCo DmLT DmF},num2cell([0.321 0.321 0.123 0.198 0.151 0.151 0.046 0.046]));
A = subs(A,{DVA MxDVA},{3.33 6.66});
A = subs(A,{CRHP CRG fMT fA},num2cell([0.25 0.25 0.25 0.6]));
A = subs(A,{DQC DQF},num2cell([36 11]*10^3));
A = subs(A,{KvHP KvG},num2cell([0.8 2]*8.7841e-06));
A = subs(A,{TauHP TauG TauMT TauVA},num2cell([0.1 0.1 1 1]*10^0));
A = subs(A,{R},num2cell(1.5*10^5));
A = subs(A,{eS},num2cell(0.6));
A = subs(A,{s s0 k cp},num2cell([1000 200 900 1000]));
A = subs(A,{VGC VR VC VMT},num2cell([19.2 133 100 0.1]*10^-3));
A = subs(A,{DdDp1 DdDp2 DdDpR DdDpC DdDpMT},num2cell([3 3 4 2.7 2.4]*10^-5)); % TODO
A = subs(A,{DdDh1 DdDh2 DdDhR DdDhC DdDhMT},num2cell([-2 -4 -1.5 -0.4 -0.3]*10^-3)); % TODO
A = subs(A,{DTDp1 DTDp2},num2cell([7.2 5.6]*10^-6)); % TODO
A = subs(A,{DTDh1 DTDh2},num2cell([4.7 55]*10^-4)); % TODO
A = subs(A,{DhDpMT},num2cell(0.01));
A = double(A);
% ----------------------- POSTPROCESSING ----------------------------------
% Normalizing with maximum deviations
% x = [DVA; p1; h1; d1; T1; p2; h2; d2; T2; DmHP; pR; hR; dR; DmG; pC; hC; dC; dMT; DmMT];
cond(A)
DVBound = 2; 
pBound = 20*10^5;
hBound = 50*10^3;
dBound = 20;
TBound = 10;
DmBound = 0.2;
T = diag(1./[DVBound; pBound; hBound; dBound; TBound; pBound; hBound; dBound; TBound; DmBound;...
    pBound; hBound; dBound; DmBound; pBound; hBound; dBound; dBound; DmBound]);
% z = Tx
% Dz = TA/Tx
A = T*A/T;
B = T*B;
G = T*G;
cond(A)
expA = sign(A).*exp(abs(A));
expA(A==0) = 0;
% -------------------------- PLOTTING -------------------------------------
figure(1)
imagesc(expA,[-10 10])
colormap('jet')
colorbar
xlabel(char(x))
ylabel(['Diff' char(x)])
title('Pieceswise exponential of normalized system A, bounded by -10...10')
% Modal analysis
% Vz = x
% Dz = V\AVz
% V\AV = D
% AV = VD
[V,eigA] = eig(A);
V*eigA/V-A
invV = inv(V);
% Separating real and imaginary values
invV = real([invV(1:9,:); real(invV(10,:)); imag(invV(10,:));...
    invV(12:14,:); real(invV(15,:)); imag(invV(15,:));...
    real(invV(17,:)); imag(invV(17,:)); invV(19,:)]);
figure(2)
subplot(131)
imagexpeigA = sign(imag(eigA)).*exp(abs(imag(eigA)));
imagexpeigA(imag(eigA)==0) = 0;
imagesc(imagexpeigA,[-10^0.5 10^0.5])
colormap('jet')
colorbar
xlabel('z')
ylabel('Dz')
title('imagexpeigA')
subplot(132)
realexpeigA = sign(real(eigA)).*exp(abs(real(eigA)));
realexpeigA(real(eigA)==0) = 0;
imagesc(realexpeigA,[-10 10])
colormap('jet')
colorbar
xlabel('z')
ylabel('Dz')
title('realexpeigA')
subplot(133)
imagesc(invV,[-10^1 10^1])
colormap('jet')
colorbar
xlabel(char(x))
ylabel('z')
title('inv(V) reordered')