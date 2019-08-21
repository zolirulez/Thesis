% Script for plotting enthalpy or temperature states (discretized to 10
% cells) in gas cooler

clearvars
addpath('matfiles')
load faultignore_ta0 
Y = h1_sim.signals.values;

% Property conversions
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
        try
            dc(it,it2-4) = CoolProp.PropsSI('D','P',Y(it,1),'H',Y(it,it2),'CO2');
        catch
            dc(it,it2-4) = CoolProp.PropsSI('D','P',Y(it,1)+1e4,'H',Y(it,it2),'CO2');
        end
    end
    hcw(it,1) = Y(it,end-9:end)*dc(it,:)'/sum(dc(it,:));
    Tcw(it,1) = CoolProp.PropsSI('T','P',Y(it,1),'H',hcw(it,1),'CO2');
end

start = 1000;
t = h1_sim.time(start:end,:);

% Time evolution of enthalpies
load faultignore
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

% Input values
U = uy_sim.signals.values(:,1:12); 
Y = uy_sim.signals.values(:,13:end); 
CRA = U(:,1);
CRIT = U(:,4);

% Delay
delay = 150;

% Time evolution of temperatures
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

% Cross covariance functions
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

% Parameter identification for dynamical model, up to order 3
Y = hcw(start:end);
X = [ones(length(Y),1) THR(start:end) TA0(start:end) CRIT(start:end) CRA(start:end)]; %  
iddata0 = iddata(Y-mean(Y),X-mean(X),1);
tfest(iddata0,0)
tfest(iddata0,1)
tfest(iddata0,2)
tfest(iddata0,3)
save('TGC_faultignore_ta0','Tcw','hcw')