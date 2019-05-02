function X = f(fFunction,XFunction,UFunction,UValue,YFunction,YValue)


for it = 1:length(UFunction)
    fFunction = subs(fFunction,{UFunction(it)},UValue(it));
end
for it = 1:length(YFunction)
    fFunction = subs(fFunction,{YFunction(it)},YValue(it));
end
for it = 1:length(XFunction)
    solve(fFunction(it),XFunction(it))
end