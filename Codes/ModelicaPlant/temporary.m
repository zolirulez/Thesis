figure, 
asdf = NaN(length(U),1);
for it = 1:length(U)-1
    asdf(it) = CoolProp.PropsSI('T','P',Y(it+1,1),'H',Y(it+1,2),'CO2');
end
hold on, plot(asdf), hold off