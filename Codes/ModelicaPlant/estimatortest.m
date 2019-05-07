clearvars
format long
% Initializing FMIKit and adding paths
addpath('C:\Users\u375749\Documents\Thesis\Codes\Linearization')
addpath('C:\Users\u375749\Documents\Thesis\Codes\MatlabPlant')
% 'Running modelcreation_simplified.m
modelcreation_simplified
% Running substitution_simplified.m
[A,B,C,D] = constantsave(A,B,C,D);
substitution_simplified
% Running kfinit_simulink.m
kfinit_simulink

% Experiment for observation
gFunction = ginit(y);
UFunction = u;
YFunction = y;
% Experiment for actuation
XFunction = x;
load uy_sim_flat

feedback = 1; % feedback of estimator
U = uy_sim.signals.values(:,1:nu);
Y = uy_sim.signals.values(:,nu+1:nu+ny); %TODO
uy = [zeros(nu+ny,1); reshape(system.A,nx*nx,1); reshape(system.B,nx*nu,1); reshape(system.C,ny*nx,1); reshape(system.D,ny*nu,1)];
ABCD = reshape([system.A system.B; system.C system.D],(nx+ny)*(nx+nu),1);
Xs = initial.xs;
finish = 10000;
start = 1001;
% Constraints (note: feedback)
w = Xs(1)*1000*1.25/(5000 + Xs(1)*6000)*2;
TBP = CoolProp.PropsSI('T','P',Y(start,1),'H',Y(start,2),'CO2');
TA1c = 1/w*TBP+(w-1)/w*U(start,10);
DmQc = 36e3/(440e3-U(start,9))+11e3/(440e3-U(start,9));
dRc = CoolProp.PropsSI('D','P',Y(start,3),'H',Y(start,4),'CO2');
Y(start,6:8) = [TA1c; DmQc; dRc];
record = NaN(nx+nx+nx+ny,finish-start+1);
for it = start:finish
    if ~rem(it,100)
        disp(['Iteration ' num2str(it)])
    end
    % ----- Estimation based on delayed parameters
    xf = kf.markovPredictor(uy(1:nu),uy(nu+1:nu+ny),reshape(uy(nu+ny+1:nu+ny+nx*nx),nx,nx),reshape(uy(nu+ny+nx*nx+1:nu+ny+nx*nx+nx*nu),nx,nu),reshape(uy(nu+ny+nx*nx+nx*nu+1:nu+ny+nx*nx+nx*nu+ny*nx),ny,nx), reshape(uy(nu+ny+nx*nx+nx*nu+ny*nx+1:nu+ny+nx*nx+nx*nu+ny*nx+ny*nu),ny,nu));
    Xs = kf.xf + Xs;
    Xs = max(0,Xs);
    UXY = [U(it,:)'; Xs; Y(it,1:ny)'];
    % ----- Set point for new iteration
    ABCD = LTVsystemDescription(UXY(1:nu), UXY(nu+1:nu+nx), UXY(nu+nx+1:nu+nx+ny),feedback);
    A = ABCD(1:nx*nx);
    B = ABCD(nx*nx+1:nx*nx+nx*nu);
    C = ABCD(nx*nx+nx*nu+1:nx*nx+nx*nu+ny*nx);
    D = ABCD(nx*nx+nx*nu+ny*nx+1:nx*nx+nx*nu+ny*nx+ny*nu);
    % Setting the new deviation points
    u = U(it+1,:)' - U(it,:)';
    % Constraints (note: feedback)
    w = Xs(1)*1000*1.25/(5000 + Xs(1)*6000)*2;
    TBP = CoolProp.PropsSI('T','P',Y(it+1,1),'H',Y(it+1,2),'CO2');
    TA1c = 1/w*TBP+(w-1)/w*U(it+1,10);
    DmQc = 36e3/(440e3-U(it+1,9))+11e3/(440e3-U(it+1,9));
    dRc = CoolProp.PropsSI('D','P',Y(it+1,3),'H',Y(it+1,4),'CO2');
    Y(it+1,6:8) = [TA1c; DmQc; dRc];
    y = Y(it+1,1:ny)' - g(gFunction,XFunction,Xs,UFunction,U(it,:));
    uy = [u; y; A; B; C; D];
    % Recording
    statecorrection = kf.Kx*kf.e;
    eigP1 = eig(kf.P1);
    record(:,it-start+1) = [Xs; eigP1; statecorrection; kf.e];
end


plotting