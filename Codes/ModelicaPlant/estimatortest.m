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
load uy_sim_faulty

U = uy_sim.signals.values(:,1:nu); % TODO
Y = uy_sim.signals.values(:,nu+1:nu+ny-1); % TODO
Y = [Y(:,1:4) 260*ones(length(Y),1) Y(:,5)];%TODO
ny = size(Y,2);
uy = [zeros(nu+ny,1); reshape(system.A,nx*nx,1); reshape(system.B,nx*nu,1); reshape(system.C,ny*nx,1); reshape(system.D,ny*nu,1); reshape(noise.Q,nx*nx,1)];
ABCDQ = [reshape([system.A system.B; system.C system.D],(nx+ny)*(nx+nu),1); reshape(noise.Q,nx*nx,1)];
Xs = initial.xs;
finish = 10000;
start = 2001;
delay = 300;
% Delay of fan and compressors
U(start-1:end,1) = U(start-1-delay:end-delay,1);
U(start-1:end,4) = U(start-1-delay:end-delay,4);
U(:,12) = 0.2*ones(length(U),1);
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
dRc = CoolProp.PropsSI('D','P',Y(start,3),'H',Y(start,4),'CO2');
d1 = CoolProp.PropsSI('D','P',Y(start,1),'H',Xs(2),'CO2');
Y(start,6) = dRc;
% RLS
rlsInitial.t = [5000; 6000; 1e3; 1e4; 1e3]*[1 1e-3 10]; 
rlsInitial.P = diag(max(rlsInitial.t')')/10; 
rls = RecursiveLeastSquares;
lambda = 1;
weight = eye(3);
rls.initialize(lambda,weight,rlsInitial);
% Low pass filter for hBP
hBP = Y(start,2);
TA0 = U(start,10);
T1 = Y(start,6);
record = NaN(nx+nx+nx+ny+nw,finish-start+1);
for it = start:finish
    if ~rem(it,100)
        disp(['Iteration ' num2str(it)])
    end
    % ----- Estimation based on delayed parameters
    xf = kf.markovPredictor(uy(1:nu),uy(nu+1:nu+ny),reshape(uy(nu+ny+1:nu+ny+nx*nx),nx,nx),reshape(uy(nu+ny+nx*nx+1:nu+ny+nx*nx+nx*nu),nx,nu),reshape(uy(nu+ny+nx*nx+nx*nu+1:nu+ny+nx*nx+nx*nu+ny*nx),ny,nx), reshape(uy(nu+ny+nx*nx+nx*nu+ny*nx+1:nu+ny+nx*nx+nx*nu+ny*nx+ny*nu),ny,nu), reshape(uy(nu+ny+nx*nx+nx*nu+ny*nx+ny*nu+1:nu+ny+nx*nx+nx*nu+ny*nx+ny*nu+nx*nx),nx,nx));
    Xs = kf.x1 + Xs;
    Xs(3) = d1; % TODO
    Xs(7) = dRc; % TODO
    % ----- Parameter estimation
%     DV = DV*(1-Ts/TauValues(2)) + Ts/TauValues(2)*U(it-1,1)*DVValues(2);
%     DmV = U(it,3)*KvValues*sqrt(U(it,6)*(Y(it,1) - Y(it,3)));
%     DQ = DmV*(U(it,11) - hBP);
%     DT = TBP - TA0;
%     sigma = TauValues(4)/Ts*(DQ/DT - sigma)+ sigma;
%     res = sigma - [1 DV]*W;
%     schur = 1 + [1 DV]*P*[1 DV]';
%     K = P*[1 DV]'/schur ;
%     W = [sigma*0.2; sigma*0.8/DV] + K*res; % [5000; 6000]
%     P = P - K*schur*K' + eye(2)/(diag(flip(W))+1e2*eye(2));
    % Param est again
    THR = CoolProp.PropsSI('T','H',U(it,11),'P',Y(it,1),'CO2');
    DmG = U(it-1,4)*VValues(3)*U(it-1,7);
    CRA = U(it-1,1);
    DV = CRA*DVValues(2);
    DmV = U(it,3)*KvValues*sqrt(U(it,6)*(Y(it,1) - Y(it,3)));
    DQ = DmV*(U(it,11) - hBP);
    DT = max(1,TBP - TA0);
    phi = [1; CRA; TA0; DmG; THR];
    out = [DQ DT Xs(8)];
    rls.regression(phi,out);
%     phiDQ'*thDQ/(phiDT'*thDT) %TODO
    W = [DQ-phi(2)*rls.t(2)'; rls.t(2)]/DT;
%     W = (DQ - )/(DT - )rls.t(1:2,1));%rls.t(1:2,1)/(phi'*rls.t(:,2));
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
%     w = DV*sigmaValues(3)*dValues(5)/(W(1) + DV*W(2));
    TBP = CoolProp.PropsSI('T','P',Y(it+1,1),'H',Y(it+1,2),'CO2');
%     DmMTc = 36e3/(445e3-U(it+1,9))+11e3/(445e3-U(it+1,9));
    dRc = CoolProp.PropsSI('D','P',Y(it+1,3),'H',Y(it+1,4),'CO2');
    d1 = CoolProp.PropsSI('D','P',Y(it+1,1),'H',Xs(2),'CO2');
    Y(it+1,5) = [dRc];
    % Low pass filters for noisy measurements
    TauNoise = 2;
    hBP = hBP*(1-Ts/TauNoise) + Ts/TauNoise*Y(it+1,2);
    TA0 = TA0*(1-Ts/TauNoise) + Ts/TauNoise*U(it+1,10); 
    T1 = T1*(1-Ts/TauNoise) + Ts/TauNoise*Y(it+1,6); 
    Y(it+1,2) = hBP;
    U(it+1,10) = TA0;
    Y(it+1,6) = T1 - (TBP - TA0);
    y = Y(it+1,1:ny)' - g(gFunction,XFunction,Xs,UFunction,U(it,:));
    kf.x1 = zeros(nx,1);
    uy = [u; y; A; B; C; D; Q];
    % Recording
    statecorrection = kf.Kx*kf.e;
    eigP1 = eig(kf.P1);
    record(:,it-start+1) = [Xs; eigP1; statecorrection; kf.e; W];
end


plotting