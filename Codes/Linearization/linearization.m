clearvars
close all
syms KvV KvG                           % Kv value of valve
syms TauV TauVA TauBP TauBC TauR TauTA TauIT  % Time constants
syms TauFake                            % Fake time constant for derivative conversions
syms R                                  % Hydraulic resistance
syms Vc VR VIT                          % Volumes of cells, reseiver and piston
syms eS                                 % Isentropic efficiency
syms p1 p2 pR                           % Pressures
syms h1 h2 hL hG hR hMT hHR hC hF       % Enthalpies
syms d1 d2 dR dG dA dBP                 % Densities
syms CRV CRA BPR BP                     % Capacity or division ratios
syms fIT                                % Frequencies
syms DmV DmG DmL Dm21                   % Mass flow rates
syms DVA MxDVA                          % Volume flow rate of air and its maximum
syms s s0 k                             % Convection parameters
syms cp                                 % Specific heat capacity of air
syms TA1 TA2 T1 T2 TA0                  % Temperatures
syms DQC DQF                            % Heat flow rates
% Partial derivatives
syms DpDd1 DpDd2 DpDdR                  % From density to pressure
syms DpDh1 DpDh2 DpDhR                  % From enthalpy to pressure
syms DTDd1 DTDd2                        % From density to temperature
syms DTDh1 DTDh2                        % From enthalpy to temperature
% Delta values
syms delta_hpIT                         % Isentropic curve delta ratio
syms delta_hHR delta_h2 delta_hJ        % Enthalpy drops
% ----------------------- STATIC EQUATIONS --------------------------------
% Boundary condition
DmC = DQC/(hC-hL);
DmF = DQF/(hF-hL);
DmMT = DmL;
% Compressor
hIT = hG + delta_hpIT*(p1 - pR)/eS;
% Joint
DmHR = DmMT + DmG;
hJ = (hMT*DmMT+hIT*DmG)/DmHR;
% Heat recovery
hHR = hJ - delta_hHR + delta_hJ;
% Gas cooler, heat transfer
s = s0 + k*DVA;
w = dA*DVA*cp/s;
DQ1 = (TA2 - T1)*s;
DQ2 = (TA1 - T2)*s;
% Gas cooler, fluid side
Joining1 = [-1 d1; 1 -DpDh1];
Joining2 = [-1 d2; 1 -DpDh2];
DmBP = DmV*BP;
Dm1 = DmHR - DmBP;
Dm2 = DmV*(1-BP);
DPsi1 = Dm1*hHR-Dm21*h1 + DQ1;
DPsi2 = Dm21*h1-Dm2*h2 + DQ2;
% ByPass valve
hBP = (Dm2*(h2-delta_h2) + DmBP*hHR)/DmV;
% High pressure valve
% Receiver
JointR = [-1 dR; 1 -DpDhR];
DPsiR = DmV*hBP - DmL*hL - DmG*hG;
% Receiver Valve
% Fan
% ---------------------- DYNAMIC EQUATIONS --------------------------------
% Boundary condition
DDmL = 1/TauBC*(-DmL + DmC + DmF);
% Heat Recovery (no differentials)
% Gas cooler
DDm21 = 1/TauR*(-Dm21 + 1/R*sqrt(d1*(p1-p2)));
Dd1 = 1/Vc*(Dm1-Dm21); % unctrb integrator, any of density dampings reduces the number by one
Dd2 = 1/Vc*(Dm21-DmV*(1-BP));
Dph1 = 1/TauFake*(-[p1; h1]+Joining1\[DPsi1/Vc; DpDd1*Dd1]); % unobsv integrator
Dph2 = 1/TauFake*(-[p2; h2]+Joining2\[DPsi2/Vc; DpDd2*Dd2]); % unctrb integrator
DT1 = 1/TauFake*(-T1+[DTDd1 DTDh1]*[Dd1; Dph1(2)]); % unctrb integrator
DT2 = 1/TauFake*(-T2+[DTDd2 DTDh2]*[Dd1; Dph2(2)]); % unctrb integrator
DTA1 = 1/TauTA*(-TA1 + 1/(w+1)*T2 + 1/(1/w+1)*TA0);
DTA2 = 1/TauTA*(-TA2 + 1/(w+1)*T1 + 1/(1/w+1)*TA1);
% ByPass Valve
DBP = 1/TauBP*(-BP + BPR);
% High pressure valve
DDmV = 1/TauV*(-DmV + CRV*KvV*sqrt(dBP*(p2 - pR)));
% Receiver
DdR = 1/VR*(DmV-DmL-DmG);
DphR = 1/TauFake*(-[pR; hR]+JointR\[DPsiR/VR; DpDdR*DdR]); % unctrb integrator
% IT Compressor
DDmG = 1/TauIT*(-DmG + dG*VIT*fIT);
% Fan (A refers to air)
DDVA = 1/TauVA*(-DVA + MxDVA*CRA);
% Disturbances
Ddelta_hJ = 1/TauFake*(-delta_hJ); % unobsv, unctrb integrator
Ddelta_h2 = 0; % unobsv, unctrb integrator
% ------------------------- MEASUREMENTS ----------------------------------
p2m = p2;
T2m = T2;
hBPm = hBP;
pRm = pR;
hRm = hR;
hHRm = hHR;
% ------------------------ LINEARIZATION ----------------------------------
% Augmentation
Dx = [DDVA; Dph1; Dd1; DT1; DTA2; DDm21; Dph2; Dd2; DT2; DTA1; DBP; DDmV; DphR; DdR; DDmG; DDmL; Ddelta_hJ; Ddelta_h2];
x = [DVA; p1; h1; d1; T1; TA2; Dm21; p2; h2; d2; T2; TA1; BP; DmV; pR; hR; dR; DmG; DmL; delta_hJ; delta_h2];
y = [p2m; T2m; hBPm; pRm; hRm; hHRm];
u = [CRA; BPR; CRV; fIT];
BC = [TA0; dA; delta_hHR; hMT; DQC; DQF];       % Decided in simulator, assumed to be known
dDer = [dBP; dG; hG; hL; hC; hF; delta_hpIT];   % Derived from estimations / measurements
dDis = [delta_hJ; delta_h2];                    % Unknown input to be estimated
d = [BC; dDis];                                 % dDer is not considered, as it is changing slowly or slightly
% Jacobians
A = jacobian(Dx,x);
B = jacobian(Dx,u);
GBC = jacobian(Dx,BC);
GdDer = jacobian(Dx,dDer);
GdDis = jacobian(Dx,dDis);
Gd = jacobian(Dx,d);
C = jacobian(y,x);
D = jacobian(y,u);
EBC = jacobian(y,BC);
EdDer = jacobian(y,dDer);
EdDis = jacobian(y,dDis);
Ed = jacobian(y,d);
% ------------------------ SUBSTITUTIONS ----------------------------------
mx4sub = struct('A',A,'B',B,'C',C,'Gd',Gd);
fields = fieldnames(mx4sub);
substitution;
A = double(mx4sub.A);
B = double(mx4sub.B);
Gd = double(mx4sub.Gd);
C = double(mx4sub.C);
% ------------------------- NORMALIZING -----------------------------------
% Normalizing with maximum required deviations
% x = [DVA; p1; h1; d1; T1; TA2; Dm21; p2; h2; d2; T2; TA1; BP; DmV; pR; hR; dR; DmG; DmL; delta_hJ; delta_h2];
condA = cond(A)
disp('The condition of the matrix indicates the relative sensitivies within the system')
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
    DmBound; DmBound;...
    hBound; hBound]);
% z = Tx
% Dz = TA/Tx
A = T*A/T;
B = T*B;
GBC = T*GBC;
GdDer = T*GdDer;
GdDis = T*GdDis;
Gd = T*Gd;
C = C/T;
condA = cond(A)
disp('The condition of the matrix decreased by normalizing')
expA = sign(A).*exp(abs(A));
expA(A==0) = 0;
% ------------------------ MODAL ANALYSIS ---------------------------------
% Modal analysis
% Wt A = E Wt
% Wt A / Wt = E
% Dz = E z
% Dz = Wt A / Wt z
% Wt \ Dz = A / Wt z
% z = Wt x
[~,E,W] = eig(A);
Wt = W';
norm(Wt\E*Wt-A,Inf)
disp('The smaller this norm, the more reliable the modal matrix is')
% Separating real and imaginary values
Wt = real([Wt(1,:); real(Wt(2,:)); imag(Wt(2,:));...
    Wt(4:8,:); real(Wt(9,:)); imag(Wt(9,:));...
    Wt(11,:); real(Wt(12,:)); imag(Wt(12,:)); Wt(14:end,:)]);
realexpeigA = sign(real(E)).*exp(abs(real(E)));
realexpeigA(real(E)==0) = 0;
imagexpeigA = sign(imag(E)).*exp(abs(imag(E)));
imagexpeigA(imag(E)==0) = 0;
% ----------------------- CONTROLLABILITY ---------------------------------
[AC,BC,CC,TC,KC] = ctrbf(A,B,C,1e-10);
% norm(AC-TC*A/TC,Inf)
% disp('The smaller this norm, the more reliable the controllability staircase form is')
rankCTRB = sum(KC)
disp(['The system is not fully controllable, there are ' num2str(length(AC)-rankCTRB)...
    ' poles in the uncontrollable subspace'])
EnC = round(real(eig(AC(1:(length(AC)-rankCTRB),1:(length(AC)-rankCTRB)))),10);
disp(['There are ' num2str(sum(EnC==0))...
    ' integrators in the uncontrollable subspace, and ' num2str(sum(EnC<0))...
    ' a stable pole'])
expAC = sign(AC).*exp(abs(AC));
expAC(round(AC,10)==0) = 0;
GdC = TC*Gd;
% ------------------------ OBSERVABILITY ----------------------------------
[AO,BO,CO,TO,KO] = obsvf(A,B,C,1e-10);
% norm(AO-TO*A/TO,Inf)
% disp('The smaller this norm, the more reliable the observability staircase form is')
rankOBSV = sum(KO)
disp(['The system is not fully observable, there are ' num2str(length(AO)-rankOBSV)...
    ' poles in the unobservable subspace'])
EnO = round(real(eig(AO(1:(length(AO)-rankOBSV),1:(length(AO)-rankOBSV)))),10);
disp(['There are ' num2str(sum(EnO==0))...
    ' integrators in the unobservable subspace, and ' num2str(sum(EnO<0))...
    ' a stable pole'])
expAO = sign(AO).*exp(abs(AO));
expAO(round(AO,10)==0) = 0;
GdO = TO*Gd;
% -------------------------- PLOTTING -------------------------------------
% System matrix
figure(1)
imagesc(expA,[-10 10])
colormap('jet')
colorbar
xlabel(char(x))
ylabel(['Diff' char(x)])
title('Pieceswise exponential of normalized system A, bounded by -10...10')
% Modal analysis
figure(2)
subplot(141)
imagesc(diag(imagexpeigA),[-10^0.5 10^0.5])
colormap('jet')
colorbar
xlabel('z')
ylabel('Dz')
title('imagexpeigA')
subplot(142)
imagesc(diag(realexpeigA),[-10 10])
colormap('jet')
colorbar
xlabel('z')
ylabel('Dz')
title('realexpeigA')
subplot(122)
imagesc(Wt,[-1 1])
colormap('jet')
colorbar
xlabel(char(x))
ylabel('z')
title('Wt, reordered')
% Controllability analysis
figure(3)
subplot(221)
imagesc(expAC,[-10^0.5 10^0.5])
colormap('jet')
colorbar
xlabel('z')
ylabel('Dz')
title('expAC')
subplot(223)
imagesc(CC,[-10^0 10^0])
colormap('jet')
colorbar
xlabel('z')
ylabel('y')
title('CC')
subplot(243)
imagesc(BC,[-10^-5 10^-5])
colormap('jet')
colorbar
xlabel(char(u))
ylabel('z')
title('BC')
subplot(244)
imagesc(GdC,[-10^-5 10^-5])
colormap('jet')
colorbar
xlabel(char(d))
ylabel('z')
title('GdC')
% Observability analysis
figure(4)
subplot(221)
imagesc(expAO,[-10^0.5 10^0.5])
colormap('jet')
colorbar
xlabel('z')
ylabel('Dz')
title('expAO')
subplot(223)
imagesc(CO,[-10^1 10^1])
colormap('jet')
colorbar
xlabel('z')
ylabel('y')
title('CO')
subplot(243)
imagesc(BO,[-10^-5 10^-5])
colormap('jet')
colorbar
xlabel(char(u))
ylabel('z')
title('BO')
subplot(244)
imagesc(GdO,[-10^-5 10^-5])
colormap('jet')
colorbar
xlabel(char(d))
ylabel('z')
title('GdO')