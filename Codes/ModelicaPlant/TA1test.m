record = [];
for it = 500:10000
    TBP = CoolProp.PropsSI('T','P',Y(it,1),'H',Y(start,2),'CO2');
    DT = TBP - U(it,10);
    TA1 = 1/w*TBP+(w-1)/w*U(it+1,10);
    record = [record [DT; TA1]];
end
figure(1)
subplot(211)
plot(record(1,:))
subplot(212)
plot(record(2,:)-273.15)