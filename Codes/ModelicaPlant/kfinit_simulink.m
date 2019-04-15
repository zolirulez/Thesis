nz = 4;
systemCont = ss(A,B,C,D);
Ts = 1;
systemDisc = c2d(systemCont,Ts);
system.A = systemDisc.A;
system.B = systemDisc.B;
system.C = systemDisc.C;
% x = [DVA; p1; h1; d1; T1; TA2; Dm21; p2; h2; d2; T2; TA1; BP; DmV; pR; hR; dR; DmG; DmL; delta_hJ; delta_h2];
Cz = zeros(nz,nx);
Cz(1,8) = 1; Cz(2,9) = 1; Cz(3,15) = 1; Cz(4,16) = 1;
system.Cz = Cz;
system.G = eye(nx); % This is another G matrix, connected to only stochastic content
initial.x = zeros(nx,1);
initial.P = 0.1*eye(nx);
kfType = 'timeinvariant';
kf = KalmanFilter;
horizon = 1;
kf.initialize(system,noise,initial,kfType,horizon)