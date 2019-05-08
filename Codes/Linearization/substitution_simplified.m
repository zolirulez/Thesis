% ------------------------ SUBSTITUTIONS ----------------------------------
DmValues = [0.321 0.123 0.198 0.198]; % DmV DmG DmL DmMT
pValues = [85 38]*10^5; % p1 pR
hValues = [(525+300)/2 ...
    CoolProp.PropsSI('H','P',pValues(2),'Q',0,'CO2')/1000 ...
    CoolProp.PropsSI('H','P',pValues(2),'Q',1,'CO2')/1000 ...
    280 ...
    (525*DmValues(4)+460*DmValues(2))/DmValues(1)]*10^3; % h1 hL hG hR hHR
hBPValue = 298*10^3;
delta_hValues = [0 hValues(1)-hBPValue]; % delta_hHR delta_h
dValues = ... % d1 dR dG dBP dA
    [CoolProp.PropsSI('D','P',pValues(1),'H',hValues(1),'CO2') ...
    CoolProp.PropsSI('D','P',pValues(2),'H',hValues(4),'CO2') ...
    CoolProp.PropsSI('D','P',pValues(2),'H',hValues(3),'CO2') ...
    CoolProp.PropsSI('D','P',pValues(1),'H',hBPValue,'CO2') ...
    1.25];
DVValues = [2 18711/3600]; % DVA MxDVA
fValues = 48; % MxfIT
KvValues = 0.8*3e-5;% 8.7841e-06;% KvV
TauValues = [1 1 10 1 0.1]; % TauV TauVA TauQ TauTA TauIT Taup
eValues = 0.65; % eS
sigmaValues = [5000 6000 1000]; % s0 k cp 
VValues = [19.2 133 0.05]*10^-3; % Vc VR VG
CRValues = [0 DmValues(1)/KvValues(1)/sqrt(dValues(1)*(pValues(1)-pValues(2)))...
    DVValues(1)/DVValues(2) DmValues(2)/dValues(3)/fValues(1)/VValues(3)]; % BP CRV CRA CRIT
wValue = 0.1;%dValues(5)*DVValues(1)*sigmaValues(3)/(sigmaValues(1)+sigmaValues(2)*DVValues(1));
TA0Value = 30+273.15;
TA1Value = 1/wValue*CoolProp.PropsSI('T','P',pValues(1),'H',hBPValue(1),'CO2')+(wValue-1)/wValue*TA0Value;
T1Value = TA1Value + 3;
TValues = [TA1Value T1Value TA0Value]; % TA1 T1 TA0
delta_phValues = [17 8]; % delta_ph1 delta_phR
delta_pdValues = [0.46 0.50]*10^4;% delta_pd1 delta_pdR
delta_ThValues = [7.0]*1e-4; % delta_Th1
delta_TdValues = [0.15]; % delta_Td1
delta_hpValues = ... % delta_hpIT
    (CoolProp.PropsSI('H','P',pValues(1),'S',...
    CoolProp.PropsSI('S','P',pValues(2),'H',hValues(3),'CO2'),'CO2')-hValues(3))/...
    (pValues(1)-pValues(2));

for it = 1:numel(fields)
    mx4sub.(fields{it}) = subs(mx4sub.(fields{it}),{DmV DmG DmL DmQ},num2cell(DmValues));
    mx4sub.(fields{it}) = subs(mx4sub.(fields{it}),{p1 pR},num2cell(pValues));
    mx4sub.(fields{it}) = subs(mx4sub.(fields{it}),{h1 hL hG hR hHR},num2cell(hValues));
    mx4sub.(fields{it}) = subs(mx4sub.(fields{it}),{delta_hHR delta_h},num2cell(delta_hValues));
    mx4sub.(fields{it}) = subs(mx4sub.(fields{it}),{d1 dR dG dBP dA},num2cell(dValues));
    mx4sub.(fields{it}) = subs(mx4sub.(fields{it}),{DVA MxDVA},num2cell(DVValues));
    mx4sub.(fields{it}) = subs(mx4sub.(fields{it}),{BP CRV CRA CRIT},num2cell(CRValues));
    mx4sub.(fields{it}) = subs(mx4sub.(fields{it}),{MxfIT},num2cell(fValues));
    mx4sub.(fields{it}) = subs(mx4sub.(fields{it}),{KvV},num2cell(KvValues));
    mx4sub.(fields{it}) = subs(mx4sub.(fields{it}),{TauV TauVA TauQ TauTA Taup},num2cell(TauValues));
    mx4sub.(fields{it}) = subs(mx4sub.(fields{it}),{eS},num2cell(eValues));
    mx4sub.(fields{it}) = subs(mx4sub.(fields{it}),{s0 k cp},num2cell(sigmaValues));
    mx4sub.(fields{it}) = subs(mx4sub.(fields{it}),{Vc VR VG},num2cell(VValues));
    mx4sub.(fields{it}) = subs(mx4sub.(fields{it}),{TA1 T1 TA0},num2cell(TValues));
    mx4sub.(fields{it}) = subs(mx4sub.(fields{it}),{delta_ph1 delta_phR},num2cell(delta_phValues));
    mx4sub.(fields{it}) = subs(mx4sub.(fields{it}),{delta_pd1 delta_pdR},num2cell(delta_pdValues));
    mx4sub.(fields{it}) = subs(mx4sub.(fields{it}),{delta_Th1},num2cell(delta_ThValues));
    mx4sub.(fields{it}) = subs(mx4sub.(fields{it}),{delta_Td1},num2cell(delta_TdValues));
    mx4sub.(fields{it}) = subs(mx4sub.(fields{it}),{delta_hpIT},num2cell(delta_hpValues));
end
A = double(mx4sub.A);
B = double(mx4sub.B);
G = double(mx4sub.G);
C = double(mx4sub.C);
D = double(mx4sub.D);