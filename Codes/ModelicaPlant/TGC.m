clearvars
load h1_sim_faulty
start = 1800;
Y = h1_sim.signals.values;
t = h1_sim.time(start:end,:);

P = [];
THR = [];
TBP = [];
TA0 = [];
Tcw = [];
Tc = NaN(length(Y),10);
dc = NaN(length(Y),10);
hcw = NaN(length(Y),1);
for it = 1:length(Y)
    P = [P; Y(it,1)];
    THR = [THR; Y(it,2)];
    TBP = [TBP; Y(it,3)];
    TA0 = [TA0; Y(it,4)];
    for it2 = 5:size(Y,2)
        Tc(it,it2-4) = CoolProp.PropsSI('T','P',Y(it,1),'H',Y(it,it2),'CO2');
        dc(it,it2-4) = CoolProp.PropsSI('D','P',Y(it,1),'H',Y(it,it2),'CO2');
    end
    Tcw = [Tcw; Tc(it,:)*dc(it,:)'/sum(dc(it,:))];
    hcw(it,1) = CoolProp.PropsSI('H','P',Y(it,1),'T',Tcw(end),'CO2');
end
%%
load uy_sim_faulty
U = uy_sim.signals.values(:,1:nu); 
Y = uy_sim.signals.values(:,nu+1:nu+ny); 
delay = 500;
THRd = THR(start-delay:end-delay);
CRA = U(start:end,1);
CRIT = U(start:end,4);
figure(5)
subplot(211)
plot(t,THRd-273,t,TBP(start:end)-273,t,TA0(start:end)-273,t,Tcw(start:end)-273)
xlabel('Time [s]')
ylabel('Temperetures [C]')
legend('T_H_R delayed','T_B_P','T_A_0','T_c weighted by density')
subplot(212)
plot(t,Tc(start:end,:)-273)
xlabel('Time [s]')
ylabel('Cell temperatures [C]')
figure(6)
[c1,lags] = xcov(Tcw(start:end),THR(start:end));
[c2,~] = xcov(Tcw(start:end),CRIT);
[c3,~] = xcov(Tcw(start:end),TA0(start:end));
[c4,~] = xcov(Tcw(start:end),CRA);
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
X = [ones(length(THRd),1) THRd TA0(start:end) CRIT CRA];
Y = Tcw(start:end);
P = (X'*X)\X'*Y;
max(abs(Y-X*P))
save('TGC')