function X = ginverse(ginv,UFunction,UValue,YFunction,YValue)
% 
% syms BP BPR
% ginv = simplify(subs(ginv,{BP},BPR));
for it = 1:length(UFunction)
    ginv = subs(ginv,{UFunction(it)},UValue(it));
end
for it = 1:length(YFunction)
    ginv = subs(ginv,{YFunction(it)},YValue(it));
end

syms DVA p1 h1 d1 TA2 Dm21 p2 h2 d2 TA1 BP DmV pR hR dR DmG delta_h2 DmQ
X = [DVA; p1; h1; d1; TA2; Dm21; p2; h2; d2; TA1; BP; DmV; pR; hR; dR; DmG; delta_h2; DmQ];
solve(ginv-YValue,X)