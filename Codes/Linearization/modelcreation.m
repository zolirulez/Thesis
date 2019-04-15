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
% ------------------------- AUGMENTATION ----------------------------------
Dx = [DDVA; Dph1; Dd1; DT1; DTA2; DDm21; Dph2; Dd2; DT2; DTA1; DBP; DDmV; DphR; DdR; DDmG; DDmL; Ddelta_hJ; Ddelta_h2];
x = [DVA; p1; h1; d1; T1; TA2; Dm21; p2; h2; d2; T2; TA1; BP; DmV; pR; hR; dR; DmG; DmL; delta_hJ; delta_h2];
y = [p2m; T2m; hBPm; pRm; hRm; hHRm];
u = [CRA; BPR; CRV; fIT];
BC = [TA0; dA; delta_hHR; hMT; DQC; DQF];       % Decided in simulator, assumed to be known
dDer = [dBP; dG; hG; hL; hC; hF; delta_hpIT];   % Derived from estimations / measurements
dDis = [delta_hJ; delta_h2];                    % Unknown input to be estimated
d = [BC; dDer; dDis];
% ------------------------ LINEARIZATION ----------------------------------                              
% Jacobians
A = jacobian(Dx,x);
B = jacobian(Dx,u);
G = jacobian(Dx,d);
C = jacobian(y,x);
D = jacobian(y,u);
E = jacobian(y,d);
mx4sub = struct('A',A,'B',B,'C',C,'G',G);
fields = fieldnames(mx4sub);