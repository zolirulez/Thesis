syms KvV KvG                           % Kv value of valve
syms TauV TauVA TauBP TauBC TauR TauTA TauIT Taup  % Time constants
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
% Conversions
syms delta_ph1 delta_ph2 delta_phR      % From enthalpy to pressure
syms delta_pd1 delta_pd2 delta_pdR      % From density to pressure
syms delta_Th2                          % From enthalpy to temperature
syms delta_Td2                          % From density to temperature
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
Dd1 = 1/Vc*(Dm1-Dm21);
Dd2 = 1/Vc*(Dm21-DmV*(1-BP));
Dp1 = 1/Taup*(-p1+delta_ph1*h1+delta_pd1*d1);
Dp2 = 1/Taup*(-p2+delta_ph2*h2+delta_pd2*d2);
Dh1 = 1/Vc*(Dm1*hHR-Dm21*h1 + DQ1);
Dh2 = 1/Vc*(Dm21*h1-Dm2*h2 + DQ2);
DTA1 = 1/TauTA*(-TA1 + 1/(w+1)*T2 + 1/(1/w+1)*TA0);
DTA2 = 1/TauTA*(-TA2 + 1/(w+1)*T1 + 1/(1/w+1)*TA1);
% ByPass Valve
DBP = 1/TauBP*(-BP + BPR);
% High pressure valve
DDmV = 1/TauV*(-DmV + CRV*KvV*sqrt(dBP*(p2 - pR)));
% Receiver
DdR = 1/VR*(DmV-DmL-DmG);
DpR = 1/Taup*(-pR+delta_phR*hR+delta_pdR*dR);
DhR = 1/VR*(DmV*hBP - DmL*hL - DmG*hG);
% IT Compressor
DDmG = 1/TauIT*(-DmG + dG*VIT*fIT);
% Fan (A refers to air)
DDVA = 1/TauVA*(-DVA + MxDVA*CRA);
% Disturbances
Ddelta_hJ = 1/TauFake*(-delta_hJ); % unobsv, unctrb integrator
Ddelta_h2 = 0; % unobsv, unctrb integrator
% ------------------------- MEASUREMENTS ----------------------------------
p2m = p2;
T2m = delta_Th2*h2+delta_Td2*d2;
hBPm = hBP;
pRm = pR;
hRm = hR;
hHRm = hHR;
% ------------------------- AUGMENTATION ----------------------------------
Dx = [DDVA; Dp1; Dh1; Dd1; DTA2; DDm21; Dp2; Dh2; Dd2; DTA1; DBP; DDmV; DpR; DhR; DdR; DDmG; DDmL; Ddelta_hJ; Ddelta_h2];
x = [DVA; p1; h1; d1; TA2; Dm21; p2; h2; d2; TA1; BP; DmV; pR; hR; dR; DmG; DmL; delta_hJ; delta_h2];
y = [p2m; T2m; hBPm; pRm; hRm; hHRm];
u = [CRA; BPR; CRV; fIT];
BC = [TA0; dA; delta_hHR; hMT; DQC; DQF];       % Decided in simulator, assumed to be known
dDer = [dBP; dG; hG; hL; hC; hF; delta_hpIT];   % Derived from estimations / measurements
dDis = [delta_hJ; delta_h2];                    % Unknown input to be estimated
d = [BC; dDer; dDis];