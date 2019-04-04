nx = length(A);
nu = size(B,2);
ny = size(C,1);
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
system.G = eye(nx);
% Noises TODO, Gd is also not included
noise.R = 0.01*eye(ny);
noise.Q = 0.01*eye(nx);
noise.S = zeros(nx,ny);
initial.x = 0.1*ones(nx,1);
initial.P = 0.1*eye(nx);
kfType = 'timeinvariant';
kf = KalmanFilter;
horizon = 100;
kf.initialize(system,noise,initial,kfType,horizon)
j = horizon;
u = zeros(nu,j);
Lr = chol(noise.R);
y = kf.C*initial.x + Lr*randn(ny,1);
Q = cell(j,1);
for it = 1:j
    Q{it} = noise.Q;
end
[xf, x1, xj, z] = kf.outputPredictor(u,y,Q,j)
[xf, x1, z] = kf.markovPredictor(u,y)
z = reshape(z,nz,j);
% Plotting
figure(1)
plot(z')
disp(['Final value should be for all-step inputs'])
kf.Cz*((eye(nx)-kf.A)\kf.B)
save('kalmanFilter','kf')