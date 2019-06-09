figure, plot(U(:,6))
dBPf = NaN(length(U),1);
for it = 1:length(U)-1
    dBPf(it) = CoolProp.PropsSI('D','P',Yf(it+1,1),'H',Yf(it+1,2),'CO2');
end
hold on, plot(dBPf), hold off