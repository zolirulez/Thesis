nz = 4;
systemCont = ss(A,B,C,D);
Ts = 1;
systemDisc = c2d(systemCont,Ts);
system.A = systemDisc.A;
system.B = systemDisc.B;
system.C = systemDisc.C;
system.D = systemDisc.D;
% x = [DVA; p1; h1; d1; TA1; BP; DmV; pR; hR; dR; DmG; delta_h2; DmQ];
% u = [CRA; BPR; CRV; CRIT; delta_hHR; dBP; dG; hG; hL; TA0; hMT];
Cz = zeros(nz,nx);
Cz(1,2) = 1; Cz(2,3) = 1; Cz(2,11) = -1; Cz(3,9) = 1; Cz(4,8) = 1;
system.Cz = Cz;
system.G = eye(nx); % This is another G matrix, connected to only stochastic content
initial.xs = [DVValues(1); pValues(1); hValues(1); dValues(1); TA1Value;...
    DmValues(1); pValues(2); hValues(4); dValues(2); DmValues(2);...
    delta_hValues(2); DmValues(4); 530e3];
initial.x = zeros(nx,1);
initial.us = [CRValues(3); CRValues(1); CRValues(2); CRValues(4); ...
    delta_hValues(1); dValues(4); dValues(3); hValues(3); hValues(2); TA0Values; hValues(5)];
initial.P = 3*noise.Q;
kfType = 'timevarying';
kf = KalmanFilter;
horizon = 1;
kf.initialize(system,noise,initial,kfType,horizon)