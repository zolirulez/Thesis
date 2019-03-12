DdDh1=CoolProp.PropsSI('d(D)/d(H)|P','P',86.5e5,'H',450e3,'CO2');
DmG=0.123;
delta_hpIT=(CoolProp.PropsSI('H','P',86.5e5,'S',...
CoolProp.PropsSI('S','P',38e5,'H',430e3,'CO2'),'CO2')-430e3)/...
(86.5e5-38e5);
DmL=0.198;
BP=0;
DmV=0.321;
Vc=19.2e-3/2;
eS=0.6;
DdDp1=CoolProp.PropsSI('d(D)/d(P)|H','P',86.5e5,'H',525e3,'CO2');
d1=CoolProp.PropsSI('D','P',86.5e5,'H',450e3,'CO2');
 -(DdDh1*DmG^2*delta_hpIT + DdDh1*DmG*DmL*delta_hpIT - BP*DdDh1*DmG*DmV*delta_hpIT)/...
     (Vc*eS*(DdDh1 + DdDp1*d1)*(DmG + DmL))
 % it should be 3.4271e+06
 (DdDp1*DmG^2*delta_hpIT + DdDp1*DmG*DmL*delta_hpIT - BP*DdDp1*DmG*DmV*delta_hpIT)/...
     (Vc*eS*(DdDh1 + DdDp1*d1)*(DmG + DmL))
 % WITHOUT pessure feedback
  -(DdDh1*DmG^2*delta_hpIT + DdDh1*DmG*DmL*delta_hpIT - BP*DdDh1*DmG*DmV*delta_hpIT)/(DdDp1*Vc*d1*eS*(DmG + DmL))
                          (DmG^2*delta_hpIT + DmG*DmL*delta_hpIT - BP*DmG*DmV*delta_hpIT)/(Vc*d1*eS*(DmG + DmL))
