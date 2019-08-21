% This script converts the measurement values, so that the faults are
% presented in temperature domain

if ~exist('fielddata')
    U = uy_sim.signals.values(:,1:nu);
    Y = uy_sim.signals.values(:,nu+1:nu+ny);
    outcorrecord = squeeze(yhat_sim.signals.values);
end
faultest = NaN(length(outcorrecord),2);
for it2 = 1:length(faultest)
    TBPf = CoolProp.PropsSI('T','P',Y(it2,1),'H',Y(it2,2),'CO2');
    try
        faultest(it2,1) = TBPf - CoolProp.PropsSI('T','P',Y(it2,1),'H',outcorrecord(2,it2),'CO2');
    catch
        faultest(it2,1) = NaN;
    end
    try
        faultest(it2,2) = TBPf - CoolProp.PropsSI('T','P',Y(it2,1),'H',outcorrecord(3,it2),'CO2');
    catch
        faultest(it2,2) = NaN;
    end
end