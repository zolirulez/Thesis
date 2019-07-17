% Script for plotting enthalpy or temperature states (discretized to 10
% cells) in gas cooler

clearvars
addpath('matfiles')
load h1_sim_faultyestctrl3 %h1_faulty_chirp % h1_sim_faulty
Y = h1_sim.signals.values;

P = [];
THR = [];
TBP = [];
TA0 = [];
Tc = NaN(length(Y),10);
dc = NaN(length(Y),10);
hcw = NaN(length(Y),1);
Tcw = NaN(length(Y),1);
hHR = NaN(length(Y),1);
for it = 1:length(Y)
    P = [P; Y(it,1)];
    THR = [THR; Y(it,2)];
    TBP = [TBP; Y(it,3)];
    TA0 = [TA0; Y(it,4)];
    hHR(it,1) = CoolProp.PropsSI('H','P',Y(it,1),'T',Y(it,2),'CO2');
    for it2 = 5:size(Y,2)
%         Tc(it,it2-4) = CoolProp.PropsSI('T','P',Y(it,1),'H',Y(it,it2),'CO2');
        dc(it,it2-4) = CoolProp.PropsSI('D','P',Y(it,1),'H',Y(it,it2),'CO2');
    end
%     Tcw = [Tcw; Tc(it,:)*dc(it,:)'/sum(dc(it,:))];
%     hcw(it,1) = CoolProp.PropsSI('H','P',Y(it,1),'T',Tcw(end),'CO2');
    hcw(it,1) = Y(it,end-9:end)*dc(it,:)'/sum(dc(it,:));
    Tcw(it,1) = CoolProp.PropsSI('T','P',Y(it,1),'H',hcw(it,1),'CO2');
end
%%
start = 1000;
t = h1_sim.time(start:end,:);

load h1_sim_faultyestctrl3 % h1_sim_faulty % h1_faulty_chirp
Y = h1_sim.signals.values;
figure(3)
subplot(211)
plot(t,hcw(start:end,:))
xlabel('Time [s]')
ylabel('Enthalpies [C]')
subplot(212)
plot(t,Y(start:end,end-9:end))
xlabel('Time [s]')
ylabel('Cell enthalpies [kJ/kg]')

load uy_sim_faultyestctrl3 % uy_sim_faulty % uy_sim_faulty_chirp
U = uy_sim.signals.values(:,1:12); 
Y = uy_sim.signals.values(:,13:end); 
CRA = U(:,1);
CRIT = U(:,4);

delay = 300;

figure(5)
subplot(211)
plot(t,THR(start-delay:end-delay)-273,t,TBP(start:end)-273,t,TA0(start:end)-273,t,Tcw(start:end)-273)
xlabel('Time [s]')
ylabel('Temperetures [C]')
legend('T_H_R delayed','T_B_P','T_A_0','T_c weighted by density')
subplot(212)
plot(t,Tc(start:end,:)-273)
xlabel('Time [s]')
ylabel('Cell temperatures [C]')

% X = [ones(length(start:length(Y)),1) THR(start-delay:end-delay) TA0(start:end) CRIT(start-delay:end-delay) CRA(start-delay:end-delay)]; %  
% Y = hcw(start:end);
% P = (X'*X)\X'*Y;
% max(abs(Y-X*P))
% 
% figure(2)
% plot(1:length(Y),hcw,t,X*P)
% P

figure(6)
[c1,lags] = xcov(THR(start-delay:end-delay),hcw(start:end));
[c2,~] = xcov(TA0(start:end),hcw(start:end));
[c3,~] = xcov(CRIT(start:end),hcw(start:end));
[c4,~] = xcov(CRA(start:end),hcw(start:end));
clf
hold on
plot(lags,-c1)
plot(lags,-c2)
plot(lags,-c3)
plot(lags,-c4)
hold off
title('Cross covariances')
xlabel('Lag [s]')
grid on

%%
Y = hcw(start:end);
X = [ones(length(Y),1) THR(start:end) TA0(start:end) CRIT(start:end) CRA(start:end)]; %  
iddata0 = iddata(Y-mean(Y),X-mean(X),1);
tfest(iddata0,0)
tfest(iddata0,1)
tfest(iddata0,2)
tfest(iddata0,3)
save('TGC_faultyestctrl3','Tcw','hcw')%save('TGC_faultyestctrl')