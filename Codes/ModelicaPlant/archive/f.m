function DX = f(DXFunction,XValue,UValue)

constants = load('constants.mat');
constantfields = fieldnames(constants);
for i=1:length(constantfields)
    eval([constantfields{i} '=constants.' constantfields{i} ';']);
end

% Input values
% u = [CRA; CRV; CRIT; CRG; DmQ; dBP; dG; hG; hL; hHR; pMT; TA0];
CRAValue = UValue(1);
CRVValue = UValue(2);
CRITValue = UValue(3);
CRGValue = UValue(4);
DmQValue = UValue(5);
dBPValue = UValue(6);
dGValue = UValue(7);
hGValue = UValue(8);
hLValue = UValue(9);
hHRValue = UValue(10);
pMTValue = UValue(11);
TA0Value = UValue(12);

% States MODIFIED (TODO)
p1Value = XValue(1);
h1Value = XValue(2);
d1Value = XValue(3);
TA1Value = XValue(4);
pRValue = XValue(5);
hRValue = XValue(6);
dRValue = XValue(7);
delta_hValue = XValue(8);
% Tables
p1idx = max(1,[find(p1Value < Pbig)-1 51]);
pRidx = max(1,[find(pRValue < Pbig)-1 51]);
h1idx = max(1,[find(h1Value < Hbig)-1 51]);
hRidx = max(1,[find(hRValue < Hbig)-1 51]);
delta_ph1Value = paramvectorP(1,p1idx(1),h1idx(1));
delta_phRValue = paramvectorP(1,pRidx(1),hRidx(1));
delta_pd1Value = paramvectorP(2,p1idx(1),h1idx(1));
delta_pdRValue = paramvectorP(2,pRidx(1),hRidx(1));
delta_Th1Value = paramvectorT(1,p1idx(1),h1idx(1));
delta_Td1Value = paramvectorT(2,p1idx(1),h1idx(1));


value = [strrep('CRA,CRG,CRIT,CRV,DmQ,TA0,TA1,d1,dBP,dG,dR,delta_h,delta_Td1,delta_Th1,delta_pd1,delta_ph1,delta_pdR,delta_phR,h1,hG,hHR,hL,hR,p1,pMT,pR',',','Value,') 'Value'];
DX = eval(['DXFunction(' value ');']);
% for it = 1:length(UFunction)
%     DXFunction = subs(DXFunction,{UFunction(it)},UValue(it));
% end
% for it = 1:length(YFunction)
%     DXFunction = subs(DXFunction,{YFunction(it)},YValue(it));
% end
% for it = 1:length(XFunction)
%     solve(DXFunction(it),XFunction(it))
% end