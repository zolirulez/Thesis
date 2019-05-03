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
feedback = 1; % feedback of estimator
U = uy_sim.signals.values(:,1:nu);
Y = uy_sim.signals.values(:,nu+1:nu+ny); %TODO
uy = [zeros(nu+ny,1); reshape(system.A,nx*nx,1); reshape(system.B,nx*nu,1); reshape(system.C,ny*nx,1); reshape(system.D,ny*nu,1)];
ABCD = reshape([system.A system.B; system.C system.D],(nx+ny)*(nx+nu),1);
Xs = initial.xs;
finish = 10000;
start = 501;
record = zeros(nx+nx+nx+ny,finish-start+1);
% TEMPORARILY TODO
Y(:,8) = Xs(6);
Y(:,9) = Xs(18) + Xs(16);
for it = start:finish
    if ~rem(it,100)
        disp(['Iteration ' num2str(it)])
    end
    % Estimation
    xf = kf.markovPredictor(uy(1:nu),uy(nu+1:nu+ny),reshape(uy(nu+ny+1:nu+ny+nx*nx),nx,nx),reshape(uy(nu+ny+nx*nx+1:nu+ny+nx*nx+nx*nu),nx,nu),reshape(uy(nu+ny+nx*nx+nx*nu+1:nu+ny+nx*nx+nx*nu+ny*nx),ny,nx), reshape(uy(nu+ny+nx*nx+nx*nu+ny*nx+1:nu+ny+nx*nx+nx*nu+ny*nx+ny*nu),ny,nu));
    Xs = xf + Xs;
    % X [DVA; p1; h1; d1; TA2; Dm21; p2; h2; d2; TA1; BP; DmV; pR; hR; dR; DmG; delta_h2; DmQ]
    % Y [p2m; hBPm; pRm; hRm; hHRm; TA1m; DmQm]
    % U [CRA; BPR; CRV; CRIT; delta_hHR; dBP; dG; hG; hL; TA0; hMT]
%     Xs([2,13:14]) = 0.1*Y(it,[1,3:4])' + 0.9*Xs([2,13:14]);
    UXY = [U(it,:)'; Xs; Y(it,:)'];
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
        TA1c = 1/w*TBP+(w-1)/w*U(it+1,10);
        DmQc = 36e3/(440e3-U(it+1,9))+11e3/(440e3-U(it+1,9));
        Dm21c = DmQc + Xs(16);
        DmVc = Dm21c;
        Y(it+1,6:9) = [TA1c; DmQc; DmVc; Dm21c];
        y = Y(it+1,:)' - Y(it,:)';
    end
    uy = [u; y; A; B; C; D];
    for itnorm = 1:nx
        normKx(itnorm,1) = norm(kf.Kx(itnorm,:));
        normP1(itnorm,1) = norm(kf.P1(itnorm,:));
    end
    record(:,it-start+1) = [Xs; normP1; normKx; kf.e];
end
t = start:Ts:finish;
figure(1)
clf
subplot(321)
hold on
plot(t,record(2,:)/1e5)
plot(t,record(7,:)/1e5)
plot(t,record(13,:)/1e5)
hold off
xlabel('Time [s]')
ylabel('Pressure [bar]')
subplot(322)
hold on
plot(t,record(3,:)/1000)
plot(t,record(8,:)/1000)
plot(t,record(14,:)/1000)
plot(t,record(17,:)/1000)
hold off
xlabel('Time [s]')
ylabel('Enthalpy [J/kg]')
subplot(323)
hold on
plot(t,record(5,:)-273.15)
plot(t,record(10,:)-273.15)
hold off
xlabel('Time [s]')
ylabel('Temperature [C]')
subplot(324)
hold on
plot(t,record(6,:))
plot(t,record(12,:))
plot(t,record(16,:))
plot(t,record(18,:))
hold off
xlabel('Time [s]')
ylabel('Mass flow rate [kg/s]')
subplot(325)
hold on
plot(t,record(4,:))
plot(t,record(9,:))
plot(t,record(15,:))
hold off
xlabel('Time [s]')
ylabel('Density [kg/m^3]')
subplot(326)
hold on
plot(t,record(1,:))
plot(t,record(11,:))
hold off

figure(2)
subplot(311)
plot(t,record(nx+1:nx+nx,:)./sqrt(diag(noise.Q)))
xlabel('Time [s]')
ylabel('Relative norm of P_1')
legend('DVA','p1','h1','d1','TA2','Dm21','p2','h2','d2','TA1','BP','DmV','pR','hR','dR','DmG','delta_h2','DmQ');
subplot(312)
plot(t,record(nx+nx+1:nx+nx+nx,:)./sqrt(diag(noise.Q)))
xlabel('Time [s]')
ylabel('Relative norm of K_x')
legend('DVA','p1','h1','d1','TA2','Dm21','p2','h2','d2','TA1','BP','DmV','pR','hR','dR','DmG','delta_h2','DmQ');
subplot(313)
plot(t,record(nx+nx+nx+1:nx+nx+nx+ny,:)./sqrt(diag(noise.R)))
legend('p_2','h_B_P','p_R','h_R','h_H_R','T_A_1','Dm_Q','Dm_V','Dm_2_1')
xlabel('Time [s]')
ylabel('Relative innovation')