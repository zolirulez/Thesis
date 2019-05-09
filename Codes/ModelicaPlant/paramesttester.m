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
Y = uy_sim.signals.values(:,nu+1:nu+ny);
W = [sigmaValues(1); sigmaValues(2)];
sigma = 74300/3;
nw = length(W);
P = diag(flip(W))/10;
start = 1000;
DV = U(start,1)*DVValues(2);
W = [5000; 6000];
Winit = W;
TAU = TauValues(2);
w = DV*sigmaValues(3)*dValues(5)/(W(1) + DV*W(2));
TBP = CoolProp.PropsSI('T','P',Y(start,1),'H',Y(start,2),'CO2');
TA1c = 1/w*TBP+(w-1)/w*U(start,10);
DmQc = 36e3/(440e3-U(start,9))+11e3/(440e3-U(start,9));
dRc = CoolProp.PropsSI('D','P',Y(start,3),'H',Y(start,4),'CO2');
Y(start,5:7) = [TA1c; DmQc; dRc];
delay = 300;
record = [];
record2 = [];
for it = start:10000
    % ----- Parameter estimation
    DV = DV*(1-Ts/TAU) + Ts/TAU*U(it-1-delay,1)*DVValues(2);
    DmV = U(it,3)*KvValues*sqrt(U(it,6)*(Y(it,1) - Y(it,3)));
    DQ = DmV*(U(it,11) - Y(it,2));
    DT = TBP - U(it,10);
    sigma = TauValues(4)/Ts*(DQ/DT - sigma)+ sigma;
    res = sigma - [1 DV]*W;
    schur = 1 + [1 DV]*P*[1 DV]';
    K = P*[1 DV]'/schur ;
    W = W + K*res;
%     W = Winit + K*res;
    P = P - K*schur*K' + eye(2)/(diag(flip(W))+1e2*eye(2));
    
    w = DV*sigmaValues(3)*dValues(5)/(W(1) + DV*W(2));
    TBP = CoolProp.PropsSI('T','P',Y(it+1,1),'H',Y(it+1,2),'CO2');
    TA1c = 1/w*TBP+(w-1)/w*U(it+1,10);
    DmQc = 36e3/(440e3-U(it+1,9))+11e3/(440e3-U(it+1,9));
    dRc = CoolProp.PropsSI('D','P',Y(it+1,3),'H',Y(it+1,4),'CO2');
    Y(it+1,5:7) = [TA1c; DmQc; dRc];
    record = [record [W; DQ; sigma; w*1e5]];
    record2 = [record2 [DV*1e4; sigma]];
end
figure(1)
subplot(211)
plot(record')
legend('s0','k','DQ','sigma','w*1e5')
ylim([-1e4 8e4])
grid
subplot(212)
plot(record2')
legend('DV*1e4','sigma')
figure(2)
N = length(record2);
[y,lags] = xcov(record2(1,:),record2(2,:));
plot(lags,y)
grid
