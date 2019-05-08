syms KvV                                % Kv value of valve
syms TauV TauVA TauQ TauTA TauIT Taup   % Time constants
syms TauFake                            % Fake time constant for derivative conversions
syms R                                  % Hydraulic resistance
syms Vc VR VG                           % Volumes of cells, reseiver and piston
syms eS                                 % Isentropic efficiency
syms p1 pR                              % Pressures
syms h1 hL hG hR hHR                    % Enthalpies
syms d1 dR dG dA dBP                    % Densities
syms CRV CRA CRIT BP                    % Capacity or division ratios
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
DmG = dG*VG*CRIT*MxfIT;
DmHR = DmMT + DmG;
% hJ = (hMT*DmMT+hIT*DmG)/DmHR;
% Heat recovery
% hHR = hJ - delta_hHR;
% Gas cooler, heat transfer
DVA = MxDVA*CRA;
s = s0 + k*DVA;
w = dA*DVA*cp/s;
T1 = delta_Th1*h1 + delta_Td1*d1;
DQ1 = (TA1 - T1)*s;
% High pressure valve
DmV = CRV*KvV*sqrt(dBP*(p1 - pR));
% Gas cooler, fluid side
DmBP = DmV*BP;
Dm1 = DmHR - DmBP;
Dm2 = DmV*(1-BP);
% ByPass valve
hBP = (Dm2*(h1-delta_h) + DmBP*hHR)/DmV;
% Receiver
% Receiver Valve
% Fan
% ---------------------- DYNAMIC EQUATIONS --------------------------------
% Heat Recovery
% Gas cooler
Dd1 = 1/Vc*(Dm1-Dm2);
Dp1 = 1/Taup*(-p1+delta_ph1*h1+delta_pd1*d1);
Dh1 = 1/(d1*Vc)*(Dm1*hHR-Dm1*hBP + DQ1);
DTA1 = 1/TauTA*(-TA1 + 1/(w+1)*T1 + 1/(1/w+1)*TA0);
% High pressure valve
% Receiver
DdR = 1/VR*(DmV-DmL-DmG);
DpR = 1/Taup*(-pR+delta_phR*hR+delta_pdR*dR);
DhR = 1/(dR*VR)*(DmV*hBP - DmL*hL - DmG*hG);
% IT Compressor
% Fan (A refers to air)
% Disturbances
Ddelta_h = 0;
DDmQ = 1/TauQ*(-DmQ + DmV - DmG);
% ------------------------- MEASUREMENTS ----------------------------------
p1m = p1;
hBPm = hBP;
pRm = pR;
hRm = hR;
TA1m = TA1;
DmQm = DmQ; % Boundary condition
dRm = dR;
% ------------------------- AUGMENTATION ----------------------------------
Dx = [Dp1; Dh1; Dd1; DTA1; DpR; DhR; DdR; Ddelta_h; DDmQ];
x = [p1; h1; d1; TA1; pR; hR; dR; delta_h; DmQ];
y = [p1m; hBPm; pRm; hRm; TA1m; DmQ; dR]; 
u = [CRA; BP; CRV; CRIT; delta_hHR; dBP; dG; hG; hL; TA0; hHR];
d = [delta_h; DmQ]; % Unknown input to be estimated
% Dimensions
nx = length(x);
nu = length(u);
ny = length(y);
% Noises
DVBound = 0.5;
pBound = 5*10^5;
hBound = 20*10^3;
dBound = 20;
TBound = 5;
DmBound = 0.1;
noise.R = diag([pBound; hBound; pBound; hBound; TBound; DmBound; dBound]);
Qcont = diag([pBound; hBound; dBound; TBound; pBound; hBound; dBound;...
    hBound; DmBound])/1e6; % TODO
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