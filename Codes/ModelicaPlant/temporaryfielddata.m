TBPtrue = NaN(length(Y),1);
TBPest = NaN(length(Y),1);
TBPestDQ = NaN(length(Y),1);
for it=7800:length(Y)
    TBPest(it) = CoolProp.PropsSI('T','P',Y(it,1),'H',outcorrecord(2,it-start),'CO2')-273;
    TBPestDQ(it) = CoolProp.PropsSI('T','P',Y(it,1),'H',Yf(it,2),'CO2')-273;
    TBPtrue(it) = CoolProp.PropsSI('T','P',Y(it,1),'H',Y(it,2),'CO2')+5-273;
end
figure,
hold on, plot(TBPtrue)
hold on, plot(TBPest)
hold on, plot(TBPestDQ)