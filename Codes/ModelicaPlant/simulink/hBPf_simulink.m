function hBPf = hBPf_simulink(DQcor,CRV,hHR,dBP,pGC,pR)

dBPf = dBP;
CRV = max(0.01,CRV);
for it = 1:3
    DmV = CRV*0.8*3e-5*sqrt(dBPf*(pGC - pR));
    hBPf = hHR-DQcor/DmV;
    try
        dBPf = CoolProp.PropsSI('D','P',pGC,'H',hBPf,'CO2');
    catch
        dBPf = CoolProp.PropsSI('D','P',pGC+1e4,'H',hBPf,'CO2');
    end
end

