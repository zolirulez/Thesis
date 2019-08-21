function q = x2q(x,p)
% Function to calculate quality from liquid height and pressure
dV = CoolProp.PropsSI('D','P',p,'Q',1,'CO2');
dL = CoolProp.PropsSI('D','P',p,'Q',0,'CO2');
q = dV*(1-x)/(dV*(1-x)+dL*x);