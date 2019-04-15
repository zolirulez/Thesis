% ------------------------ SUBSTITUTIONS ----------------------------------
pValues = [86.5 85 38]*10^5; % p1 p2 pR
hValues = [450 325 212 430 280 525 525]*10^3; % h1 h2 hL hG hR hMT hHR
hBPValue = 298*10^3;
delta_hValues = [0 0 hValues(2)-hBPValue]; % delta_hHR delta_hJ delta_h2
dValues = ... % d1 d2 dR dG dBP dA
    [CoolProp.PropsSI('D','P',pValues(1),'H',hValues(1),'CO2') ...
    CoolProp.PropsSI('D','P',pValues(2),'H',hValues(2),'CO2') ...
    CoolProp.PropsSI('D','P',pValues(3),'H',hValues(5),'CO2') ...
    CoolProp.PropsSI('D','P',pValues(3),'H',hValues(4),'CO2') ...
    CoolProp.PropsSI('D','P',pValues(2),'H',hBPValue,'CO2') ...
    1.2];
DmValues = [0.321 0.123 0.198 0.321 0.198]; % DmV DmG DmL Dm21 DmMT
DVValues = [3.33 6.66]; % DVA MxDVA
CRValues = [0 0.25 0.6 0.22]; % BP CRV CRA CRIT
fValues = 48; % MxfIT
KvValues = [0.8 2]*8.7841e-06;% KvV KvG
TauValues = [1 5 1 10 1 10 5 0.1]; % TauV TauVA TauBP TauQ TauR TauTA TauIT Taup
RValues = 2.5*10^4; % R
eValues = 0.6; % eS
sigmaValues = [5000/2 6000/2 1000]; % s0 k cp
VValues = [19.2/2 133 0.05]*10^-3; % Vc VR VIT
wValues = dValues(6)*DVValues(1)*sigmaValues(3)/(sigmaValues(1)+sigmaValues(2)*DVValues(1));
T1Values = CoolProp.PropsSI('T','P',pValues(1),'H',hValues(1),'CO2');
T2Values = CoolProp.PropsSI('T','P',pValues(2),'H',hValues(2),'CO2');
TA0Values = 30+273.15;
TA1Values = 1/(wValues+1)*T2Values+wValues/(wValues+1)*TA0Values;
TA2Values = 1/(wValues+1)*T1Values+wValues/(wValues+1)*TA1Values;
TValues = [TA2Values TA1Values T1Values T2Values TA0Values]; % TA2 TA1 T1 T2 TA0
delta_pdValues = [1.25 0.46 0.50]*10^4;% delta_pd1 delta_pd2 delta_pdR
delta_phValues = [13 17 8]; % delta_ph1 delta_ph2 delta_phR
delta_TdValues = [0.16 0.15]; % delta_Td1 delta_Td2
delta_ThValues = [6.6 7.0]*1e-4; % delta_Th1 delta_Th2
delta_hpValues = ... % delta_hpIT
    (CoolProp.PropsSI('H','P',pValues(1),'S',...
    CoolProp.PropsSI('S','P',pValues(3),'H',hValues(4),'CO2'),'CO2')-hValues(4))/...
    (pValues(1)-pValues(3));

for it = 1:numel(fields)
    mx4sub.(fields{it}) = subs(mx4sub.(fields{it}),{p1 p2 pR},num2cell(pValues));
    mx4sub.(fields{it}) = subs(mx4sub.(fields{it}),{h1 h2 hL hG hR hMT hHR},num2cell(hValues));
    mx4sub.(fields{it}) = subs(mx4sub.(fields{it}),{delta_hHR delta_hJ delta_h2},num2cell(delta_hValues));
    mx4sub.(fields{it}) = subs(mx4sub.(fields{it}),{d1 d2 dR dG dBP dA},num2cell(dValues));
    mx4sub.(fields{it}) = subs(mx4sub.(fields{it}),{DmV DmG DmL Dm21 DmQ},num2cell(DmValues));
    mx4sub.(fields{it}) = subs(mx4sub.(fields{it}),{DVA MxDVA},num2cell(DVValues));
    mx4sub.(fields{it}) = subs(mx4sub.(fields{it}),{BP CRV CRA CRIT},num2cell(CRValues));
    mx4sub.(fields{it}) = subs(mx4sub.(fields{it}),{MxfIT},num2cell(fValues));
    mx4sub.(fields{it}) = subs(mx4sub.(fields{it}),{KvV KvG},num2cell(KvValues));
    mx4sub.(fields{it}) = subs(mx4sub.(fields{it}),{TauV TauVA TauBP TauQ TauR TauTA TauIT Taup},num2cell(TauValues));
    mx4sub.(fields{it}) = subs(mx4sub.(fields{it}),{R},num2cell(RValues));
    mx4sub.(fields{it}) = subs(mx4sub.(fields{it}),{eS},num2cell(eValues));
    mx4sub.(fields{it}) = subs(mx4sub.(fields{it}),{s0 k cp},num2cell(sigmaValues));
    mx4sub.(fields{it}) = subs(mx4sub.(fields{it}),{Vc VR VIT},num2cell(VValues));
    mx4sub.(fields{it}) = subs(mx4sub.(fields{it}),{TA2 TA1 T1 T2 TA0},num2cell(TValues));
    mx4sub.(fields{it}) = subs(mx4sub.(fields{it}),{delta_pd1 delta_pd2 delta_pdR},num2cell(delta_pdValues));
    mx4sub.(fields{it}) = subs(mx4sub.(fields{it}),{delta_ph1 delta_ph2 delta_phR},num2cell(delta_phValues));
    mx4sub.(fields{it}) = subs(mx4sub.(fields{it}),{delta_Td1 delta_Td2},num2cell(delta_TdValues));
    mx4sub.(fields{it}) = subs(mx4sub.(fields{it}),{delta_Th1 delta_Th2},num2cell(delta_ThValues));
    mx4sub.(fields{it}) = subs(mx4sub.(fields{it}),{delta_hpIT},num2cell(delta_hpValues));
end
A = double(mx4sub.A);
B = double(mx4sub.B);
G = double(mx4sub.G);
C = double(mx4sub.C);
D = double(mx4sub.D);