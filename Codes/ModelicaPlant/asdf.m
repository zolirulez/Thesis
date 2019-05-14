TBP = [];
for it = 1:length(Y)
TBP = [TBP CoolProp.PropsSI('T','P',Y(it,1),'H',Y(it,2),'CO2')];
end
plot(TBP-273)
hold on
plot(U(:,10)-273)
hold off