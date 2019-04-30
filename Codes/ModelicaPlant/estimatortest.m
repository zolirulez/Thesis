clearvars
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

U = uy_sim.signals.values(:,1:nu);
Y = uy_sim.signals.values(:,nu+1:nu+ny);
u = U - initial.us';
y = Y - (kf.C*initial.xs + kf.D*initial.us)';
uy = [zeros(nu+ny,1); reshape(system.A,nx*nx,1); reshape(system.B,nx*nu,1); reshape(system.C,ny*nx,1); reshape(system.D,ny*nu,1)];
for it = 1000:2000
    xf = kf.markovPredictor(uy(1:nu),uy(nu+1:nu+ny),reshape(uy(nu+ny+1:nu+ny+nx*nx),nx,nx),reshape(uy(nu+ny+nx*nx+1:nu+ny+nx*nx+nx*nu),nx,nu),reshape(uy(nu+ny+nx*nx+nx*nu+1:nu+ny+nx*nx+nx*nu+ny*nx),ny,nx), reshape(uy(nu+ny+nx*nx+nx*nu+ny*nx+1:nu+ny+nx*nx+nx*nu+ny*nx+ny*nu),ny,nu));
    Xf = xf + initial.xs;
    ux = [u(it,:)'+initial.us; xf+initial.xs];
    ABCD = LTVsystemDescription(ux(1:nu),ux(nu+1:nu+nx));
    A = ABCD(1:nx*nx);
    B = ABCD(nx*nx+1:nx*nx+nx*nu);
    C = ABCD(nx*nx+nx*nu+1:nx*nx+nx*nu+ny*nx);
    D = ABCD(nx*nx+nx*nu+ny*nx+1:nx*nx+nx*nu+ny*nx+ny*nu);
    uy = [u(it,:)'; y(it,:)'; A; B; C; D];
end