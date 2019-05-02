syms KvV KvG                           % Kv value of valve
syms TauV TauVA TauBP TauQ TauR TauTA TauIT Taup  % Time constants
syms TauFake                            % Fake time constant for derivative conversions
syms R                                  % Hydraulic resistance
syms Vc VR VG                          % Volumes of cells, reseiver and piston
syms eS                                 % Isentropic efficiency
syms p1 p2 pR                           % Pressures
syms h1 h2 hL hG hR hMT hMTd hHR        % Enthalpies
syms d1 d2 dR dG dA dBP                 % Densities
syms CRV CRA CRIT BPR BP                % Capacity or division ratios
syms MxfIT                              % Frequencies
syms DmV DmG DmL Dm21 DmQ               % Mass flow rates
syms DVA MxDVA                          % Volume flow rate of air and its maximum
syms s s0 k                             % Convection parameters
syms cp                                 % Specific heat capacity of air
syms TA1 TA2 TA0                        % Air temperatures
% Conversions
syms delta_ph1 delta_ph2 delta_phR      % From enthalpy to pressure
syms delta_pd1 delta_pd2 delta_pdR      % From density to pressure
syms delta_Th1 delta_Th2                % From enthalpy to temperature
syms delta_Td1 delta_Td2                % From density to temperature
% Delta values
syms delta_hpIT                         % Isentropic curve delta ratio
syms delta_hHR delta_h2                 % Enthalpy drops
% ----------------------- STATIC EQUATIONS --------------------------------
% Boundary condition
DmL = DmQ;
DmMT = DmQ;
% Compressor
hIT = hG + delta_hpIT*(p1 - pR)/eS;
% Joint
DmHR = DmMT + DmG;
hJ = (hMT*DmMT+hIT*DmG)/DmHR;
% Heat recovery
hHR = hJ - delta_hHR;
% Gas cooler, heat transfer
s = s0 + k*DVA;
w = dA*DVA*cp/s;
T1 = delta_Th1*h1 + delta_Td1*d1;
T2 = delta_Th2*h2 + delta_Td2*d2;
DQ1 = (TA2 - T1)*s;
DQ2 = (TA1 - T2)*s;
% Gas cooler, fluid side
DmBP = DmV*BP;
Dm1 = DmHR - DmBP;
Dm2 = DmV*(1-BP);
% ByPass valve
hBP = (Dm2*(h2-delta_h2) + DmBP*hHR)/DmV;
% High pressure valve
% Receiver
% Receiver Valve
% Fan
% ---------------------- DYNAMIC EQUATIONS --------------------------------
% Heat Recovery (no differentials)
% Gas cooler
DDm21 = 1/TauR*(-Dm21 + 1/R*sqrt(d1*(p1-p2)));
Dd1 = 1/Vc*(Dm1-Dm21);
Dd2 = 1/Vc*(Dm21-Dm2);
Dp1 = 1/Taup*(-p1+delta_ph1*h1+delta_pd1*d1);
Dp2 = 1/Taup*(-p2+delta_ph2*h2+delta_pd2*d2);
Dh1 = 1/(d1*Vc)*(Dm1*hHR-Dm21*h1 + DQ1);
Dh2 = 1/(d2*Vc)*(Dm21*h1-Dm2*h2 + DQ2);
DTA1 = 1/TauTA*(-TA1 + 1/(w+1)*T2 + 1/(1/w+1)*TA0);
DTA2 = 1/TauTA*(-TA2 + 1/(w+1)*T1 + 1/(1/w+1)*TA1);
% ByPass Valve
DBP = 1/TauBP*(-BP + BPR);
% High pressure valve
DDmV = 1/TauV*(-DmV + CRV*KvV*sqrt(dBP*(p2 - pR)));
% Receiver
DdR = 1/VR*(DmV-DmL-DmG);
DpR = 1/Taup*(-pR+delta_phR*hR+delta_pdR*dR);
DhR = 1/(dR*VR)*(DmV*hBP - DmL*hL - DmG*hG);
% IT Compressor
DDmG = 1/TauIT*(-DmG + dG*VG*CRIT*MxfIT);
% Fan (A refers to air)
DDVA = 1/TauVA*(-DVA + MxDVA*CRA);
% Disturbances
Ddelta_h2 = 0;
DDmQ = 0;
% ------------------------- MEASUREMENTS ----------------------------------
p2m = p2;
hBPm = hBP;
pRm = pR;
hRm = hR;
hHRm = hHR;
TA1m = TA1;
DmQm = DmQ; % Boundary condition
DmVm = DmV; % Should add up
Dm21m = Dm21;
h1m = h1;
% ------------------------- AUGMENTATION ----------------------------------
Dx = [DDVA; Dp1; Dh1; Dd1; DTA2; DDm21; Dp2; Dh2; Dd2; DTA1; DBP; DDmV; DpR; DhR; DdR; DDmG; Ddelta_h2; DDmQ];
x = [DVA; p1; h1; d1; TA2; Dm21; p2; h2; d2; TA1; BP; DmV; pR; hR; dR; DmG; delta_h2; DmQ];
y = [p2m; hBPm; pRm; hRm; hHRm; TA1m; DmQm; DmVm; Dm21m; h1m]; 
u = [CRA; BPR; CRV; CRIT; delta_hHR; dBP; dG; hG; hL; TA0; hMT];
d = [delta_h2; DmQ]; % Unknown input to be estimated
% Dimensions
nx = length(x);
nu = length(u);
ny = length(y);
% Noises
DVBound = 0.005;
pBound = 5*10^5;
hBound = 20*10^3;
BPBound = 0;
dBound = 20;
TBound = 5;
DmBound = 0.005;
noise.R = diag([pBound*0; hBound*0; pBound*0; hBound*0; hBound*0;...
    TBound; DmBound; DmBound; DmBound; hBound*1e3]);
noise.Q = diag([DVBound*0; pBound; hBound; dBound; TBound; DmBound; pBound;...
    hBound; dBound; TBound; BPBound*0; DmBound; pBound; hBound; dBound;...
    DmBound*0; hBound*1e2; DmBound*1e2])*1e2; % TODO
noise.S = zeros(nx,ny);
% ------------------------ LINEARIZATION ----------------------------------                              
% Jacobians
A = jacobian(Dx,x);
B = jacobian(Dx,u);
G = jacobian(Dx,d);
C = jacobian(y,x);
D = jacobian(y,u);
E = jacobian(y,d);
mx4sub = struct('A',A,'B',B,'C',C,'D',D,'G',G);
fields = fieldnames(mx4sub);