function Y = g(gFunction,XFunction,XValue,UFunction,UValue)

load deltaValues
syms delta_Td1 delta_Th1

for it = 1:length(UFunction)
    gFunction = subs(gFunction,{UFunction(it)},UValue(it));
end
for it = 1:length(XFunction)
    gFunction = subs(gFunction,{XFunction(it)},XValue(it));
end
gFunction = subs(gFunction,{delta_Td1},delta_Td1Value);
gFunction = subs(gFunction,{delta_Th1},delta_Th1Value);

Y = double(gFunction);