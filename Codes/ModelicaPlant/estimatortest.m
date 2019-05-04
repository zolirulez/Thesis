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

% Experiment for observation
UFunction = u;
YFunction = y;
% Experiment for actuation
XFunction = x;
load uy_sim

lti = 0;
interactionDamping = 0.8;
feedback = 1; % feedback of estimator
U = uy_sim.signals.values(:,1:nu);
Y = uy_sim.signals.values(:,nu+1:nu+ny); %TODO
uy = [zeros(nu+ny,1); reshape(system.A,nx*nx,1); reshape(system.B,nx*nu,1); reshape(system.C,ny*nx,1); reshape(system.D,ny*nu,1)];
ABCD = reshape([system.A system.B; system.C system.D],(nx+ny)*(nx+nu),1);
Xs = initial.xs*1;
finish = 10000;
start = 501;
record = NaN(nx+nx+nx+ny,finish-start+1);
for it = start:finish
    if ~rem(it,100)
        disp(['Iteration ' num2str(it)])
    end
    % Estimation
    xf = kf.markovPredictor(uy(1:nu),uy(nu+1:nu+ny),reshape(uy(nu+ny+1:nu+ny+nx*nx),nx,nx),reshape(uy(nu+ny+nx*nx+1:nu+ny+nx*nx+nx*nu),nx,nu),reshape(uy(nu+ny+nx*nx+nx*nu+1:nu+ny+nx*nx+nx*nu+ny*nx),ny,nx), reshape(uy(nu+ny+nx*nx+nx*nu+ny*nx+1:nu+ny+nx*nx+nx*nu+ny*nx+ny*nu),ny,nu));
    Xs = xf*interactionDamping + Xs;
    UXY = [U(it,:)'; Xs; Y(it,1:ny)'];
    % Set point
    if ~lti
        ABCD = LTVsystemDescription(UXY(1:nu), UXY(nu+1:nu+nx), UXY(nu+nx+1:nu+nx+ny),feedback);
    end
    A = ABCD(1:nx*nx);
    B = ABCD(nx*nx+1:nx*nx+nx*nu);
    C = ABCD(nx*nx+nx*nu+1:nx*nx+nx*nu+ny*nx);
    D = ABCD(nx*nx+nx*nu+ny*nx+1:nx*nx+nx*nu+ny*nx+ny*nu);
    if lti
        u = uvec(it,:)';
        y = yvec(it,:)';
    else
        u = U(it+1,:)' - U(it,:)';
        % Constraints (note: feedback)
        w = Xs(1)*1000*1.25/(5000 + Xs(1)*6000)*2;
        TBP = CoolProp.PropsSI('T','P',Y(it+1,1),'H',Y(it+1,2),'CO2');
        TA1c = -1/w*TBP+(w+1)/w*U(it+1,10);
        DmQc = 36e3/(440e3-U(it+1,9))+11e3/(440e3-U(it+1,9));
        Dm21c = Xs(12);
        Y(it+1,6:8) = [TA1c; DmQc; Dm21c];
        if it > start
            y = Y(it+1,1:ny)' - (Y(it,1:ny)' + (kf.C*kf.x1 + kf.D*u)*interactionDamping);
        else
            y = zeros(ny,1);
        end
    end
    uy = [u; y; A; B; C; D];
    statecorrection = kf.Kx*kf.e;
    eigP1 = eig(kf.P1);
    record(:,it-start+1) = [Xs; eigP1; statecorrection; kf.e];
end


plotting