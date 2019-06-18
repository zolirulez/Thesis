function Y = g(YFunction,XValue,UValue)

p1Value = XValue(1);
h1Value = XValue(2);
d1Value = XValue(3);
TA1Value = XValue(4);
pRValue = XValue(5);
hRValue = XValue(6);
dRValue = XValue(7);
delta_hValue = XValue(8);
value = [strrep('delta_h,h1,hR,p1,pR',',','Value,') 'Value'];
Y = eval(['YFunction(' value ');']);