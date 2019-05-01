clearvars
format long
% Initializing FMIKit and adding paths
FMIKit.initialize
addpath('C:\Users\u375749\Documents\Thesis\Codes\Linearization')
addpath('C:\Users\u375749\Documents\Thesis\Codes\MatlabPlant')
% 'Running modelcreation_simplified.m
modelcreation_simplified
% Running substitution_simplified.m
[A,B,C,D] = constantsave(A,B,C,D);
substitution_simplified
% Running kfinit_simulink.m
kfinit_simulink

load uy_sim

feedback = 1; % feedback of estimator
U = uy_sim.signals.values(:,1:nu);
Y = uy_sim.signals.values(:,nu+1:nu+ny);
% u = U - initial.us';
% y = Y - (kf.C*initial.xs + kf.D*initial.us)';
uy = [zeros(nu+ny,1); reshape(system.A,nx*nx,1); reshape(system.B,nx*nu,1); reshape(system.C,ny*nx,1); reshape(system.D,ny*nu,1)];
ABCD = reshape([system.A system.B; system.C system.D],(nx+ny)*(nx+nu),1);
Xs = initial.xs;
finish = 2000;
start = 201;
record = zeros(nx,finish-start+1);
for it = start:finish
    if ~rem(it,100)
        disp(['Iteration ' num2str(it)])
    end
    xf = kf.markovPredictor(uy(1:nu),uy(nu+1:nu+ny),reshape(uy(nu+ny+1:nu+ny+nx*nx),nx,nx),reshape(uy(nu+ny+nx*nx+1:nu+ny+nx*nx+nx*nu),nx,nu),reshape(uy(nu+ny+nx*nx+nx*nu+1:nu+ny+nx*nx+nx*nu+ny*nx),ny,nx), reshape(uy(nu+ny+nx*nx+nx*nu+ny*nx+1:nu+ny+nx*nx+nx*nu+ny*nx+ny*nu),ny,nu));
    Xs = xf + Xs;%E = kf.C*xf + kf.D*u(it,:)' + Xf;
    UXY = [U(it,:)'; Xs; Y(it,:)'];
    ABCD = LTVsystemDescription(UXY(1:nu), UXY(nu+1:nu+nx), UXY(nu+nx+1:nu+nx+ny),feedback);
    A = ABCD(1:nx*nx);
    B = ABCD(nx*nx+1:nx*nx+nx*nu);
    C = ABCD(nx*nx+nx*nu+1:nx*nx+nx*nu+ny*nx);
    D = ABCD(nx*nx+nx*nu+ny*nx+1:nx*nx+nx*nu+ny*nx+ny*nu);
    u = U(it,:)' - U(it-1,:)'; % Is it right? TODO
    y = Y(it,:)' - Y(it-1,:)';
    uy = [u; y; A; B; C; D];
    record(:,it-start+1) = Xs;
end
t = start:Ts:finish;
figure(1)
clf
subplot(221)
hold on
plot(t,record(2,:))
plot(t,record(7,:))
plot(t,record(13,:))
hold off
subplot(222)
hold on
plot(t,record(3,:))
plot(t,record(8,:))
plot(t,record(14,:))
plot(t,record(17,:))
hold off
subplot(223)
hold on
plot(t,record(5,:)-273)
plot(t,record(10,:)-273)
hold off
subplot(224)
hold on
plot(t,record(6,:))
plot(t,record(12,:))
plot(t,record(16,:))
plot(t,record(18,:))
hold off