figure, plot(Y(:,2))
hBPtrue = NaN(length(Y),1);
for it=5000:length(Y)
    hBPtrue(it) = CoolProp.PropsSI('H','P',Y(it,1),'T',CoolProp.PropsSI('T','P',Y(it,1),'H',Y(it,2),'CO2')+5,'CO2');
end
hold on, plot(hBPtrue)
hold on, plot(start:start+length(outcorrecord)-1,outcorrecord(2,:)')
hold on, plot(Yf(:,2))