clearvars
close all
syms KvV KvG                           % Kv value of valve
syms TauV TauG TauVA TauBP TauBC TauR TauTA TauIT  % Time constants
syms R                                  % Hydraulic resistance
syms Vc VR VIT                          % Volumes of cells, reseiver and piston
syms eS                                 % Isentropic efficiency
syms p1 p2 pR                           % Pressures
syms h1 h2 hL hG hR hMT hHR hC hF       % Enthalpies
syms d1 d2 dR dG dA dBP                 % Densities
syms CRV BPR BP                         % Capacity or division ratios
syms fIT fA                             % Frequencies
syms DmV DmG DmL Dm21                   % Mass flow rates
syms DVA MxDVA                          % Volume flow rate of air and its maximum
syms s s0 k                             % Convection parameters
syms cp                                 % Specific heat capacity of air
syms TA1 TA2 T1 T2 TBP TA0              % Temperatures
syms DQC DQF                            % Heat flow rates
% Partial derivatives
syms DdDp1 DdDp2 DdDpR                  % From pressure to density
syms DdDh1 DdDh2 DdDhR                  % From enthalpy to density
syms DTDp1 DTDp2                        % From pressure to temperature
syms DTDh1 DTDh2                        % From enthalpy to temperature
% Delta values
syms delta_hpIT                         % Isentropic curve delta ratio
syms delta_hHR delta_h2                 % Enthalpy drops
% ----------------------- STATIC EQUATIONS --------------------------------
% Boundary condition
DmC = DQC/(hC-hL);
DmF = DQF/(hF-hL);
DmMT = DmL;
% Compressor
hIT = hG + delta_hpIT*(p1 - pR)/eS;
% Heat recovery
DmHR = DmMT + DmG;
hHR = (hMT*DmMT+hIT*DmG)/DmHR - delta_hHR;
% Gas cooler, heat transfer
s = s0 + k*DVA;
w = dA*DVA*cp/s;
DQ1 = (TA2 - T1)*s;
DQ2 = (TA1 - T2)*s;
% Gas cooler, fluid side
Joining1 = [-1 d1; DdDp1 DdDh1];
Joining2 = [-1 d2; DdDp2 DdDh2];
DmBP = DmV*BP;
Dm1 = DmHR - DmBP;
Dm2 = DmV*(1-BP);
DPsi1 = Dm1*hHR-Dm21*h1 + DQ1;
DPsi2 = Dm21*h1-Dm2*h2 + DQ2;
% ByPass valve
hBP = (Dm2*(h2-delta_h2) + DmBP*hHR)/DmV;
% High pressure valve
% Receiver
JointR = [-1 dR; DdDpR DdDhR];
DPsiR = DmV*hBP - DmL*hL - DmG*hG;
% Receiver Valve
% Fan
% ---------------------- DYNAMIC EQUATIONS --------------------------------
% Boundary condition
DDmL = 1/TauBC*(-DmL + DmC + DmF);
% Heat Recovery (no differentials)
% Gas cooler
DDm21 = 1/TauR*(-Dm21 + 1/R*sqrt(d1*(p1-p2)));
Dd1 = 1/Vc*(Dm1-Dm21);
Dd2 = 1/Vc*(Dm21-DmV*(1-BP));
Dph1 = Joining1\[DPsi1/Vc; Dd1];
Dph2 = Joining2\[DPsi2/Vc; Dd2];
DT1 = [DTDp1 DTDh1]*Dph1;
DT2 = [DTDp2 DTDh2]*Dph2;
DTA1 = 1/TauTA*(-TA1 + 1/(w+1)*T2 + 1/(1/w+1)*TA0);
DTA2 = 1/TauTA*(-TA2 + 1/(w+1)*T1 + 1/(1/w+1)*TA1);
% ByPass Valve
DBP = 1/TauBP*(-BP + BPR);
% High pressure valve
DDmV = 1/TauV*(-DmV + CRV*KvV*sqrt(dBP*(p2 - pR)));
% Receiver
DdR = 1/VR*(DmV-DmL-DmG);
DphR = JointR\[DPsiR/VR; DdR];
% IT Compressor
DDmG = 1/TauIT*(-DmG + dG*VIT*fIT);
% Fan (A refers to air)
DDVA = 1/TauVA*(-DVA + dA*MxDVA*fA);
% ------------------------ LINEARIZATION ----------------------------------
% Augmentation
Dx = [DDVA; Dph1; Dd1; DT1; DTA2; DDm21; Dph2; Dd2; DT2; DTA1; DBP; DDmV; DphR; DdR; DDmG; DDmL];
x = [DVA; p1; h1; d1; T1; TA2; Dm21; p2; h2; d2; T2; TA1; BP; DmV; pR; hR; dR; DmG; DmL];
u = [fA; BPR; CRV; fIT];
BC = [TA0; dA; delta_hHR; hC; hF; hMT; DQC; DQF];
d = [TBP; dBP; dG; delta_h2; hG; hL; delta_hpIT];
% Jacobians
A = jacobian(Dx,x);
B = jacobian(Dx,u);
GBC = jacobian(Dx,BC);
Gd = jacobian(Dx,d);
% Substitutions check the values
A = subs(A,{p1 p2 pR},num2cell([86.5 85 38]*10^5));
A = subs(A,{h1 h2 hL hG hR hMT hHR hC hF},num2cell([450 325 212 430 290 525 525 460 460]*10^3));
A = subs(A,{delta_hHR delta_h2},num2cell([0 325-298]*10^3));
A = subs(A,{d1 d2 dR dG dBP dA},...
    num2cell([CoolProp.PropsSI('D','P',86.5e5,'H',450e3,'CO2') ...
    CoolProp.PropsSI('D','P',85e5,'H',325e3,'CO2') ...
    CoolProp.PropsSI('D','P',38e5,'H',290e3,'CO2') ...
    CoolProp.PropsSI('D','P',38e5,'H',430e3,'CO2') ...
    CoolProp.PropsSI('D','P',85e5,'H',298e3,'CO2') ...
    1.2]));
A = subs(A,{TA2 TA1 T1 T2 TBP TA0},...
    num2cell([CoolProp.PropsSI('T','P',86.5e5,'H',450e3,'CO2')-20 ...
    CoolProp.PropsSI('T','P',85e5,'H',325e3,'CO2')-3 ...
    CoolProp.PropsSI('T','P',86.5e5,'H',450e3,'CO2') ...
    CoolProp.PropsSI('T','P',85e5,'H',325e3,'CO2') ...
    33+273.15 ...
    30+273.15]));
A = subs(A,{DmV DmG DmL Dm21},num2cell([0.321 0.123 0.198 0.321]));
A = subs(A,{DVA MxDVA},{3.33 6.66});
A = subs(A,{BP CRV},num2cell([0 0.25]));
A = subs(A,{DQC DQF},num2cell([36 11]*10^3));
A = subs(A,{KvV KvG},num2cell([0.8 2]*8.7841e-06));
A = subs(A,{TauV TauG TauVA TauBP TauBC TauR TauTA TauIT},num2cell([1 1 5 1 10 1 10 5]));
A = subs(A,{R},num2cell(2.5*10^4));
A = subs(A,{eS},num2cell(0.6));
A = subs(A,{s0 k cp},num2cell([5000/2 6000/2 1000]));
A = subs(A,{Vc VR VIT},num2cell([19.2/2 133 0.1]*10^-3));
A = subs(A,{DdDp1 DdDp2 DdDpR},num2cell([...
    CoolProp.PropsSI('d(D)/d(P)|H','P',86.5e5,'H',450e3,'CO2') ...
    CoolProp.PropsSI('d(D)/d(P)|H','P',85e5,'H',325e3,'CO2') ...
    CoolProp.PropsSI('d(D)/d(P)|H','P',38e5,'H',290e3,'CO2')]));
A = subs(A,{DdDh1 DdDh2 DdDhR},num2cell([...
    CoolProp.PropsSI('d(D)/d(H)|P','P',86.5e5,'H',450e3,'CO2') ...
    CoolProp.PropsSI('d(D)/d(H)|P','P',85e5,'H',325e3,'CO2') ...
    CoolProp.PropsSI('d(D)/d(H)|P','P',38e5,'H',290e3,'CO2')]));
A = subs(A,{DTDp1 DTDp2},num2cell([...
    CoolProp.PropsSI('d(T)/d(P)|H','P',86.5e5,'H',450e3,'CO2') ...
    CoolProp.PropsSI('d(T)/d(P)|H','P',85e5,'H',325e3,'CO2')]));
A = subs(A,{DTDh1 DTDh2},num2cell([...
    CoolProp.PropsSI('d(T)/d(H)|P','P',86.5e5,'H',450e3,'CO2') ...
    CoolProp.PropsSI('d(T)/d(H)|P','P',85e5,'H',325e3,'CO2')]));
A = subs(A,{delta_hpIT},num2cell(...
    (CoolProp.PropsSI('H','P',86.5e5,'S',...
CoolProp.PropsSI('S','P',38e5,'H',430e3,'CO2'),'CO2')-430e3)/...
(86.5e5-38e5)));
A = double(A);
% ----------------------- POSTPROCESSING ----------------------------------
% Normalizing with maximum required deviations
% x = [DVA; p1; h1; d1; T1; TA2; Dm21; p2; h2; d2; T2; TA1; BP; DmV; pR; hR; dR; DmG; DmL];
condA = cond(A)
DVBound = 1;
pBound = 5*10^5;
hBound = 20*10^3;
BPBound = 0.1;
dBound = 10;
TBound = 5;
DmBound = 0.1;
T = diag(1./[DVBound;...
    pBound; hBound; dBound; TBound; TBound;...
    DmBound; ...
    pBound; hBound; dBound; TBound; TBound;...
    BPBound; DmBound; ...
    pBound; hBound; dBound;...
    DmBound; DmBound]);
% z = Tx
% Dz = TA/Tx
A = T*A/T;
B = T*B;
GBC = T*GBC;
Gd = T*Gd;
condA = cond(A)
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