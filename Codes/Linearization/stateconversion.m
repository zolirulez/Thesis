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
                        D(h,p) = D(h-1,p);
                        T(h,p) = T(h-1,p);
                    else
                        D(h,p) = D(h,p-1);
                        T(h,p) = T(h,p-1);
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
figure(1)
surf(Hbig(2:end),Pbig(2:end),squeeze(paramvectorP(1,:,:)))
xlabel('h')
ylabel('p')
title('par_p_h')
figure(2)
surf(Hbig(2:end),Pbig(2:end),squeeze(paramvectorP(2,:,:)))
xlabel('h')
ylabel('p')
title('par_p_d')
figure(3)
surf(Hbig(2:end),Pbig(2:end),squeeze(paramvectorT(1,:,:)))
xlabel('h')
ylabel('p')
title('par_T_h')
figure(4)
surf(Hbig(2:end),Pbig(2:end),squeeze(paramvectorT(2,:,:)))
xlabel('h')
ylabel('p')
title('par_T_d')
% plot(y-X*param)