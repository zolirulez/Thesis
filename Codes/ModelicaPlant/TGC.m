clearvars
load h1_sim_chaos
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
    hcw(it,it2-4) = CoolProp.PropsSI('H','P',Y(it,1),'T',Tcw(end),'CO2');
end
%%
THRd = THR;
delay = 350;
THRd = THR(start-delay:end-delay);
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
[c,lags] = xcov(Tcw(start:end),THR(start:end));
plot(lags,c)
title('Cross covariance between THR and Tcw')
xlabel('Lag [s]')
X = [ones(length(THRd),1) THRd TA0(start:end)];
Y = Tcw(start:end);
P = (X'*X)\X'*Y;
max(abs(Y-X*P))
save('TGC')