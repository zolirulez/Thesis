%% This program plots a surface of partial differentials for a ph range

clearvars
% close all
handle = figure(1);
set(handle, 'Position',  [100, 100, 100+800, 100+200])
P = linspace(30,100,50)*1e5;
H = linspace(200,525,50)*1e3;
Diff = NaN(50);
for p = 1:length(P)
    for h = 1:length(H)
        Diff(h,p) = CoolProp.PropsSI('D(P)/D(H)|D','P',P(p),'H',H(h),'CO2');
        if abs(Diff(h,p))>500
            Diff(h,p) = NaN;
        end
%         Diff(h,p) = sign(Diff(h,p))*min([abs(Diff(h,p)) 1e7]);
    end
end
subplot(121)
surf(H/1e3,P/1e5,Diff/1e5*1e3)
xlabel('h [kJ/kg]')
ylabel('p [bar]')
title('$\frac{\partial p}{\partial h}|_{\rho}$ [bar kJ$^{-1}$kg]','Interpreter','latex')
for p = 1:length(P)
    for h = 1:length(H)
        Diff(h,p) = CoolProp.PropsSI('D(P)/D(D)|H','P',P(p),'H',H(h),'CO2');
        Diff(h,p) = sign(Diff(h,p))*min([abs(Diff(h,p)) 1e7]);
    end
end
subplot(122)
for p = 1:length(P)
    for h = 1:length(H)
        Diff(h,p) = CoolProp.PropsSI('D(P)/D(D)|H','P',P(p),'H',H(h),'CO2');
        if abs(Diff(h,p))>1e7
            Diff(h,p) = NaN;
        end
%         Diff(h,p) = sign(Diff(h,p))*min([abs(Diff(h,p)) 1e7]);
    end
end
surf(H/1e3,P/1e5,Diff/1e5)
xlabel('h [kJ/kg]')
ylabel('p [bar]')
title('$\frac{\partial p}{\partial \rho}|_{h}$ [bar kg m$^{-3}$]','Interpreter','latex')
saveas(handle,'partials.png')