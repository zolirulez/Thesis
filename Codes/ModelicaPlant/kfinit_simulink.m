nz = 4;
systemCont = ss(A,B,C,D);
Ts = 1;
systemDisc = c2d(systemCont,Ts);
system.A = systemDisc.A;
system.B = systemDisc.B;
system.C = systemDisc.C;
system.D = systemDisc.D;
% x = [DVA; p1; h1; d1; TA2; Dm21; p2; h2; d2; TA1; BP; DmV; pR; hR; dR; DmG; delta_h2; hMT; TA0; DmQ];
% u = [CRA; BPR; CRV; CRIT; delta_hHR; dBP; dG; hG; hL; TA0; hMT];
Cz = zeros(nz,nx);
Cz(1,8) = 1; Cz(2,9) = 1; Cz(3,15) = 1; Cz(4,16) = 1;
system.Cz = Cz;
system.G = eye(nx); % This is another G matrix, connected to only stochastic content
initial.xs = [DVValues(1); pValues(1); hValues(1); dValues(1); TA2Values; DmValues(4);...
    pValues(2); hValues(2); dValues(2); TA1Values; CRValues(1); DmValues(1);...
    pValues(3); hValues(5); dValues(3); DmValues(2);...
    delta_hValues(2); DmValues(5)];
initial.x = zeros(nx,1);
initial.us = [CRValues(3); CRValues(1); CRValues(2); CRValues(4); ...
    delta_hValues(1); dValues(5); dValues(4); hValues(4); hValues(3); TA0Values; hValues(6)];
initial.P = 30*noise.Q;
kfType = 'timevarying';
kf = KalmanFilter;
horizon = 1;
kf.initialize(system,noise,initial,kfType,horizon)