syms KvV KvG                            % Kv value of valve
syms TauV TauVA TauQ TauTA TauIT Taup Tauh   % Time constants
syms TauFake                            % Fake time constant for derivative conversions
syms R                                  % Hydraulic resistance
syms Vc VR VG VMT                       % Volumes of cells, reseiver and piston
syms eS                                 % Isentropic efficiency
syms p1 pR pMT                          % Pressures
syms h1 hL hG hR hHR                    % Enthalpies
syms d1 dR dG dA dBP                    % Densities
syms CRV CRA CRIT CRG                   % Capacity or division ratios
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
% hIT = hG + delta_hpIT*(p1 - pR)/eS;
% Joint
DmIT = dG*VG*CRIT*MxfIT;
DmG = CRG*KvG*sqrt(dG*(pR - pMT));
DmHR = DmMT + DmG;
% hJ = (hMT*DmMT+hIT*DmG)/DmHR;
% Heat recovery
% hHR = hJ - delta_hHR;
% Gas cooler, heat transfer
DVA = MxDVA*CRA;
s = s0 + k*DVA;
w = dA*DVA*cp/s;
T1 = delta_Th1*h1 + delta_Td1*d1;
%DQ1 = (TA1 - T1)*s;
DQ1 = (TA0 - TA1)*dA*DVA*cp;
% High pressure valve
DmV = CRV*KvV*sqrt(dBP*(p1 - pR));
% Gas cooler, fluid side
% ByPass valve
hBP = h1-delta_h;
% Receiver
% Receiver Valve
% Fan
% ---------------------- DYNAMIC EQUATIONS --------------------------------
% Heat Recovery
% Gas cooler
Dd1 = 1/Vc*(DmHR-DmV);
Dp1 = 1/Taup*(-p1+delta_ph1*h1+delta_pd1*d1);
Dh1 = 1/(d1*Vc)*(DmHR*hHR-DmV*hBP + DQ1);
DTA1 = 1/TauTA*(-TA1 + 1/(w+1)*T1 + 1/(1/w+1)*TA0);
% High pressure valve
% Receiver
DdR = 1/VR*(DmV-DmL-DmIT-DmG);
DpR = 1/Taup*(-pR+delta_phR*hR+delta_pdR*dR);
DhR = 1/(dR*VR)*(DmV*hBP - DmL*hL - DmG*hG - DmIT*hG);
% IT Compressor
% Fan (A refers to air)
% Disturbances
Ddelta_h = -1/Tauh*delta_h*0;
% ------------------------- MEASUREMENTS ----------------------------------
p1m = p1;
hBPm = hBP;
pRm = pR;
hRm = hR;
h1m = h1;
% ------------------------- AUGMENTATION ----------------------------------
Dx = [Dp1; Dh1; Dd1; DTA1; DpR; DhR; DdR; Ddelta_h];
x = [p1; h1; d1; TA1; pR; hR; dR; delta_h];
y = [p1m; hBPm; pRm; hRm; h1m]; 
u = [CRA; CRV; CRIT; CRG; DmQ; dBP; dG; hG; hL; hHR; pMT; TA0];
% Dimensions
nx = length(x);
nu = length(u);
ny = length(y);
% Noises
DVBound = 0.5;
pBound = 2*10^5;
hBound = 5*10^3;
dBound = 5;
TBound = 2;
DmBound = 0.05;
% Note: simulate with correct fillingratio noise!
noise.R = diag([pBound; hBound; pBound; hBound; hBound])*1e3;
Qcont = diag([pBound*1e2; hBound*1e2; dBound; TBound;...
    pBound*1e2; hBound*1e2; dBound;...
    hBound*1e2]); % TODO
noise.S = zeros(nx,ny);
% ------------------------ LINEARIZATION ----------------------------------                              
% Jacobians
A = jacobian(Dx,x);
B = jacobian(Dx,u);
C = jacobian(y,x);
D = jacobian(y,u);
mx4sub = struct('A',A,'B',B,'C',C,'D',D);
fields = fieldnames(mx4sub);