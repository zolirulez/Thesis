clearvars
format long
% Initializing FMIKit and adding paths
addpath('C:\Users\u375749\Documents\Thesis\Codes\Linearization')
addpath('C:\Users\u375749\Documents\Thesis\Codes\MatlabPlant')
% 'Running modelcreation_simplified.m
modelcreation_simplified
% Running substitution_simplified.m
[A,B,C,D] = constantsave(A,B,C,D,Qcont);
substitution_simplified
% Running kfinit_simulink.m
kfinit_simulink

% Experiment for observation
gFunction = ginit(y);
UFunction = u;
YFunction = y;
% Experiment for actuation
XFunction = x;
load uy_sim_sin

U = uy_sim.signals.values(:,1:nu);
Y = uy_sim.signals.values(:,nu+1:nu+ny); %TODO
uy = [zeros(nu+ny,1); reshape(system.A,nx*nx,1); reshape(system.B,nx*nu,1); reshape(system.C,ny*nx,1); reshape(system.D,ny*nu,1); reshape(noise.Q,nx*nx,1)];
ABCDQ = [reshape([system.A system.B; system.C system.D],(nx+ny)*(nx+nu),1); reshape(noise.Q,nx*nx,1)];
Xs = initial.xs;
finish = 10000;
start = 2001;
delay = 300;
% Delay of fan
U(start-1:end,1) = U(start-1-delay:end-delay,1);
% Parameter estimation
W = [sigmaValues(1); sigmaValues(2)];
sigma = 74300/3;
nw = length(W);
P = diag(W)/1e3;
DV = U(start,1)*DVValues(2);
W = [5000; 6000];
% Constraints (note: feedback)
w = DV*sigmaValues(3)*dValues(5)/(W(1) + DV*W(2));
TBP = CoolProp.PropsSI('T','P',Y(start,1),'H',Y(start,2),'CO2');
DmQc = 36e3/(440e3-U(start,9))+11e3/(440e3-U(start,9));
dRc = CoolProp.PropsSI('D','P',Y(start,3),'H',Y(start,4),'CO2');
Y(start,5:6) = [DmQc; dRc];
record = NaN(nx+nx+nx+ny+nw,finish-start+1);
for it = start:finish
    if ~rem(it,100)
        disp(['Iteration ' num2str(it)])
    end
    % ----- Estimation based on delayed parameters
    xf = kf.markovPredictor(uy(1:nu),uy(nu+1:nu+ny),reshape(uy(nu+ny+1:nu+ny+nx*nx),nx,nx),reshape(uy(nu+ny+nx*nx+1:nu+ny+nx*nx+nx*nu),nx,nu),reshape(uy(nu+ny+nx*nx+nx*nu+1:nu+ny+nx*nx+nx*nu+ny*nx),ny,nx), reshape(uy(nu+ny+nx*nx+nx*nu+ny*nx+1:nu+ny+nx*nx+nx*nu+ny*nx+ny*nu),ny,nu), reshape(uy(nu+ny+nx*nx+nx*nu+ny*nx+ny*nu+1:nu+ny+nx*nx+nx*nu+ny*nx+ny*nu+nx*nx),nx,nx));
    Xs = kf.xf + Xs;
    % ----- Parameter estimation
    DV = DV*(1-Ts/TauValues(2)) + Ts/TauValues(2)*U(it-1,1)*DVValues(2);
    DmV = U(it,3)*KvValues*sqrt(U(it,6)*(Y(it,1) - Y(it,3)));
    DQ = DmV*(U(it,11) - Y(it,2));
    DT = TBP - U(it,10);
    sigma = TauValues(4)/Ts*(DQ/DT - sigma)+ sigma;
    res = sigma - [1 DV]*W;
    schur = 1 + [1 DV]*P*[1 DV]';
    K = P*[1 DV]'/schur ;
    W = [5000; 6000] + K*res; % [5000; 6000]
    P = P - K*schur*K' + eye(2)/(diag(flip(W))+1e2*eye(2));
    % ----- Set point for new iteration
    UXYW = [U(it,:)'; Xs; Y(it,1:ny)'; W];
    ABCDQ = LTVsystemDescription(UXYW(1:nu), UXYW(nu+1:nu+nx), UXYW(nu+nx+1:nu+nx+ny), UXYW(nu+nx+ny+1:nu+nx+ny+nw));
    A = ABCDQ(1:nx*nx);
    B = ABCDQ(nx*nx+1:nx*nx+nx*nu);
    C = ABCDQ(nx*nx+nx*nu+1:nx*nx+nx*nu+ny*nx);
    D = ABCDQ(nx*nx+nx*nu+ny*nx+1:nx*nx+nx*nu+ny*nx+ny*nu);
    Q = ABCDQ(nx*nx+nx*nu+ny*nx+ny*nu+1:nx*nx+nx*nu+ny*nx+ny*nu+nx*nx);
    % Setting the new deviation points
    u = U(it+1,:)' - U(it,:)';
    % Constraints (note: feedback)
    w = DV*sigmaValues(3)*dValues(5)/(W(1) + DV*W(2));
    TBP = CoolProp.PropsSI('T','P',Y(it+1,1),'H',Y(it+1,2),'CO2');
    DmQc = 36e3/(440e3-U(it+1,9))+11e3/(440e3-U(it+1,9));
    dRc = CoolProp.PropsSI('D','P',Y(it+1,3),'H',Y(it+1,4),'CO2');
    Y(it+1,5:6) = [DmQc; dRc];
    y = Y(it+1,1:ny)' - g(gFunction,XFunction,Xs,UFunction,U(it,:));
    uy = [u; y; A; B; C; D; Q];
    % Recording
    statecorrection = kf.Kx*kf.e;
    eigP1 = eig(kf.P1);
    record(:,it-start+1) = [Xs; eigP1; statecorrection; kf.e; W];
end


plotting