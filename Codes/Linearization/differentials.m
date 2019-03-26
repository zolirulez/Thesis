clearvars
close all
P = linspace(30,100,100)*1e5;
H = linspace(200,525,100)*1e3;
Diff = NaN(100);
for p = 1:length(P)
    for h = 1:length(H)
        Diff(h,p) = CoolProp.PropsSI('D(P)/D(T)|H','P',P(p),'H',H(h),'CO2');
        Diff(h,p) = sign(Diff(h,p))*min([abs(Diff(h,p)) 1e6]);
    end
end
surf(H,P,Diff)