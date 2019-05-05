function Y = g(gFunction,XFunction,XValue,UFunction,UValue)

for it = 1:length(XFunction)
    gFunction = subs(gFunction,{XFunction(it)},XValue(it));
end
for it = 1:length(UFunction)
    gFunction = subs(gFunction,{UFunction(it)},UValue(it));
end

Y = double(gFunction);