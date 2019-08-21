clearvars
close all
simulationstring = 'faultcontrol_ta0';
TsParam = 1000;
estimatortest_simulink
TA0 = U(1001:5000,12);
TBP = NaN(4000,1);
for it = 1001:5000
    TBP(it-1000) = CoolProp.PropsSI('T','P',Y(it,1),'H',Y(it,2),'CO2');
end
h = figure(1);
plot(smooth(TA0-273.15,100),smooth(TBP-273.15,100))
[c, lags] = xcorr(detrend(TBP),detrend(TA0));
figure(2);
set(h, 'Position',  [100, 100, 100+700, 100+200])
plot(lags,c)

main_meas
TA0 = U(1:7000,12);
TBP = NaN(7000,1);
for it = 1:7000
    TBP(it) = CoolProp.PropsSI('T','P',Y(it,1),'H',Y(it,2),'CO2');
end
h = figure(1)
[c, lags] = xcorr(detrend(TBP),detrend(TA0));

[c, lags] = xcorr(detrend(TBP),detrend(TA0));
h = figure(2);
hold on
plot(lags,c,'r')
hold off
grid on
grid minor
xlim([-2000 2000])
xlabel('Lag [s]')
ylabel('Cross-covariance [K$^2$]','Interpreter','latex')
legend('Simulation data','Field data')
saveas(h,'crosscov.png')

iddata0 = iddata(detrend(TBP),detrend(TA0),1);
tfest(iddata0,1)
