nz = 4;
Ts = 1;
Lq = chol(Qcont,'lower');
[A,B,Qdisc] = c2dn(A,B,Lq,Ts);
system.A = A;
system.B = B;
system.C = C;
system.D = D;
noise.Q = Qdisc;
Cz = zeros(nz,nx);
Cz(1,2) = 1; Cz(2,3) = 1; Cz(2,11) = -1; Cz(3,7) = 1; Cz(4,8) = 1;
Cz = C;
system.Cz = Cz;
system.G = eye(nx); % This is another G matrix, connected to only stochastic content
initial.xs = [pValues(1); hValues(1); dValues(1); TA1Value;...
    pValues(2); hValues(4); dValues(2);...
    delta_hValues(2); DmValues(4)];
initial.x = zeros(nx,1);
initial.us = [CRValues(3); CRValues(1); CRValues(2); CRValues(4); ...
    delta_hValues(1); dValues(4); dValues(3); hValues(3); hValues(2);...
    TA0Value; hValues(5)];
initial.P = 3*Qcont;
kfType = 'timevarying';
kf = KalmanFilter;
horizon = 1;
kf.initialize(system,noise,initial,kfType,horizon)

% Noises
sigma2T = (1/3)^2;
sigma2p = (0.5*1e5/3)^2;
sigma2x = (0.01/3)^2;