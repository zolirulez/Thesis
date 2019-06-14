THR = NaN(length(Y),1);
for it = 1:length(Y)
    THR(it) = CoolProp.PropsSI('T','P',Y(it,1),'H',U(it,10),'CO2');
end