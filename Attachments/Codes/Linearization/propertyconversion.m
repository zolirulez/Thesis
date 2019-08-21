% This script is to produce the table of LSQ transfer parameters
% Zoltan Mark Pinter, Master Thesis, 2019

% Normal dataset
clearvars
close all
resolution = 5;
division = 50;
Pbig = linspace(30,100,division+1)*1e5;
Hbig = linspace(200,525,division+1)*1e3;
maxerrorP = 0;
maxerrorT = 0;
paramvectorP = NaN(2,division-1,division-1);
paramvectorT = NaN(2,division-1,division-1);
for it1 = 1:division
    for it2 = 1:division
%         disp(['Iteration ' num2str(it1) ' , ' num2str(it2)])
        P = linspace(Pbig(it1),Pbig(it1+1),resolution);
        H = linspace(Hbig(it2),Hbig(it2+1),resolution);
        D = NaN(resolution);
        T = NaN(resolution);
        states = [];
        for p = 1:length(P)
            for h = 1:length(H)
                try
                    D(h,p) = CoolProp.PropsSI('D','P',P(p),'H',H(h),'CO2');
                    T(h,p) = CoolProp.PropsSI('T','P',P(p),'H',H(h),'CO2');
                catch
                    disp('bug')
                    if h > 1
                    D(h,p) = CoolProp.PropsSI('D','P',P(p)+1e4,'H',H(h),'CO2');
                    T(h,p) = CoolProp.PropsSI('T','P',P(p)+1e4,'H',H(h),'CO2');
                    else
                    D(h,p) = CoolProp.PropsSI('D','P',P(p)+1e4,'H',H(h),'CO2');
                    T(h,p) = CoolProp.PropsSI('T','P',P(p)+1e4,'H',H(h),'CO2');
                    end
                end
                states = [states; P(p),T(h,p),H(h),D(h,p)];
            end
        end 
        yP = states(:,1);
        yT = states(:,2);
        X = states(:,3:end);
        paramP = (X'*X)\X'*yP;
        paramT = (X'*X)\X'*yT;
        errorP = (yP-X*paramP);
        errorT = (yT-X*paramT);
        maxerrorP = max(max(abs(errorP(:,1))),maxerrorP);
        maxerrorT = max(max(abs(errorT(:,1))),maxerrorT);
        paramvectorP(:,it1,it2) = paramP;
        paramvectorT(:,it1,it2) = paramT;
    end
end
maxerrorP
maxerrorT
% Plotting
if 1
    handle = figure(1);
    set(handle, 'Position',  [100, 100, 100+800, 100+400])
    subplot(221)
    surf(Hbig(2:end)/1e3,Pbig(2:end)/1e5,squeeze(paramvectorP(1,:,:))/1e5*1e3)
    xlabel('h [kJ/kg]')
    ylabel('p [bar]')
    title('$\partial_{ph}$ [bar kJ$^{-1}$kg]','Interpreter','latex')
    subplot(222)
    surf(Hbig(2:end)/1e3,Pbig(2:end)/1e5,squeeze(paramvectorP(2,:,:))/1e5)
    xlabel('h [kJ/kg]')
    ylabel('p [bar]')
    title('$\partial_{p\rho}$ [bar kg m$^{-3}$]','Interpreter','latex')
    subplot(223)
    surf(Hbig(2:end)/1e3,Pbig(2:end)/1e5,squeeze(paramvectorT(1,:,:))*1e3)
    xlabel('h [kJ/kg]')
    ylabel('p [bar]')
    title('$\partial_{Th}$ [K kJ$^{-1}$kg]','Interpreter','latex')
    subplot(224)
    surf(Hbig(2:end)/1e3,Pbig(2:end)/1e5,squeeze(paramvectorT(2,:,:)))
    xlabel('h [kJ/kg]')
    ylabel('p [bar]')
    title('$\partial_{T\rho}$ [K kg m$^{-3}$]','Interpreter','latex')
end
saveas(handle,'transferparam.png')
% Saving
% save('transferparam.mat','Hbig','Pbig','paramvectorP','paramvectorT')