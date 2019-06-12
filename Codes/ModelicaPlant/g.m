function Y = g(YFunction,XValue,UValue)

% load deltaValues
% syms delta_Td1 delta_Th1

% for it = 1:length(UFunction)
%     gFunction = subs(gFunction,{UFunction(it)},UValue(it));
% end
% for it = 1:length(XFunction)
%     gFunction = subs(gFunction,{XFunction(it)},XValue(it));
% end
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
% gFunction = subs(gFunction,{delta_Td1},delta_Td1Value);
% gFunction = subs(gFunction,{delta_Th1},delta_Th1Value);
% if XValue(4) - ( delta_Th1Value*XValue(2)+delta_Td1Value*XValue(3))>0
%     warning('Wrong heat flow sign')
% end

% Y = double(gFunction);