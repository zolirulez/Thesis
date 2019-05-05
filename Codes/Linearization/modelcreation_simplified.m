syms KvV KvG                           % Kv value of valve
syms TauV TauVA TauBP TauQ TauTA TauIT Taup  % Time constants
syms TauFake                            % Fake time constant for derivative conversions
syms R                                  % Hydraulic resistance
syms Vc VR VG                           % Volumes of cells, reseiver and piston
syms eS                                 % Isentropic efficiency
syms p1 pR                              % Pressures
syms h1 hL hG hR hMT hHR                % Enthalpies
syms d1 dR dG dA dBP                    % Densities
syms CRV CRA CRIT BPR BP                % Capacity or division ratios
syms MxfIT                              % Frequencies
syms DmV DmG DmL DmQ                    % Mass flow rates
syms DVA MxDVA                          % Volume flow rate of air and its maximum
syms s s0 k                             % Convection parameters
syms cp                                 % Specific heat capacity of air
syms TA1 TA0                            % Air temperatures
% Conversions
syms delta_ph1 delta_phR                % From enthalpy to pressure
syms delta_pd1 delta_pdR                % From density to pressure
syms delta_Th1                          % From enthalpy to temperature
syms delta_Td1                          % From density to temperature
% Delta values
syms delta_hpIT                         % Isentropic curve delta ratio
syms delta_hHR delta_h                  % Enthalpy drops
% ----------------------- STATIC EQUATIONS --------------------------------
% Boundary condition
DmL = DmQ;
DmMT = DmQ;
% Compressor
hIT = hG + delta_hpIT*(p1 - pR)/eS;
% Joint
DmHR = DmMT + DmG;
hJ = (hMT*DmMT+hIT*DmG)/DmHR;
% Gas cooler, heat transfer
s = s0 + k*DVA;
w = dA*DVA*cp/s;
T1 = delta_Th1*h1 + delta_Td1*d1;
DQ1 = (TA1 - T1)*s;
% Gas cooler, fluid side
DmBP = DmV*BP;
Dm1 = DmHR - DmBP;
Dm2 = DmV*(1-BP);
% ByPass valve
hBP = (Dm2*(h1-delta_h) + DmBP*hHR)/DmV;
% High pressure valve
% Receiver
% Receiver Valve
% Fan
% ---------------------- DYNAMIC EQUATIONS --------------------------------
% Heat Recovery
DhHR = 1/TauQ*(-hHR + hJ - delta_hHR);
% Gas cooler
Dd1 = 1/Vc*(Dm1-Dm2);
Dp1 = 1/Taup*(-p1+delta_ph1*h1+delta_pd1*d1);
Dh1 = 1/(d1*Vc)*(Dm1*hHR-Dm1*h1 + DQ1);
DTA1 = 1/TauTA*(-TA1 + 1/(w+1)*T1 + 1/(1/w+1)*TA0);
% ByPass Valve
DBP = 1/TauBP*(-BP + BPR);
% High pressure valve
DDmV = 1/TauV*(-DmV + CRV*KvV*sqrt(dBP*(p1 - pR)));
% Receiver
DdR = 1/VR*(DmV-DmL-DmG);
DpR = 1/Taup*(-pR+delta_phR*hR+delta_pdR*dR);
DhR = 1/(dR*VR)*(DmV*hBP - DmL*hL - DmG*hG);
% IT Compressor
DDmG = 1/TauIT*(-DmG + dG*VG*CRIT*MxfIT);
% Fan (A refers to air)
DDVA = 1/TauVA*(-DVA + MxDVA*CRA);
% Disturbances
Ddelta_h = 0;
DDmQ = 0;
% ------------------------- MEASUREMENTS ----------------------------------
p1m = p1;
hBPm = hBP;
pRm = pR;
hRm = hR;
hHRm = hHR;
TA1m = TA1;
DmQm1 = DmQ; % Boundary condition
DmQm2 = DmQ; % Mass balance
dRm = dR;
% ------------------------- AUGMENTATION ----------------------------------
Dx = [DDVA; Dp1; Dh1; Dd1; DTA1; DBP; DDmV; DpR; DhR; DdR; DDmG; Ddelta_h; DDmQ; DhHR];
x = [DVA; p1; h1; d1; TA1; BP; DmV; pR; hR; dR; DmG; delta_h; DmQ; hHR];
y = [p1m; hBPm; pRm; hRm; hHRm; TA1m; DmQm1; DmQm2; dR]; 
u = [CRA; BPR; CRV; CRIT; delta_hHR; dBP; dG; hG; hL; TA0; hMT];
d = [delta_h; DmQ]; % Unknown input to be estimated
% Dimensions
nx = length(x);
nu = length(u);
ny = length(y);
% Noises
DVBound = 0.1;
pBound = 5*10^5;
hBound = 20*10^3;
BPBound = 1e-4;
dBound = 20;
TBound = 5;
DmBound = 0.0005;
noise.R = diag([pBound*1e-4; hBound; pBound*1e-4; hBound*1e-4; hBound*1e-4;...
    TBound*1e-4; DmBound; DmBound; dBound*1e-4]);
noise.Q = diag([DVBound*1e-4; pBound; hBound; dBound; TBound; BPBound*1e-4;...
    DmBound*1e-4; pBound; hBound; dBound; DmBound*1e-4; hBound*1e2; DmBound*1e2; hBound*1e2]); % TODO
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