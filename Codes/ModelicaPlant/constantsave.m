function [A,B,C,D] = constantsave(A,B,C,D)

symvars = symvar([A B; C D]);
symvarsstr = ['[', strjoin(arrayfun(@char, symvars, 'uniform', 0),', '), ']'];
endindex = find(any([symvarsstr == ','; symvarsstr == ']']));
startindex = find(any([symvarsstr == '['; symvarsstr == ' ']));
value = cell(length(symvars),1);
for it = 1:length(symvars)
    value{it} = [symvarsstr(startindex(it)+1:endindex(it)-1) 'Value'];
end

% Table required
delta_ph1Value = 13;
delta_ph2Value = 17;
delta_phRValue = 8;
delta_pd1Value = 1.25*1e4;
delta_pd2Value = 0.46*1e4;
delta_pdRValue = 0.50*1e4;
delta_Th1Value = 6.6*1e-4;
delta_Th2Value = 7.0*1e-4;
delta_Td1Value = 0.16;
delta_Td2Value = 0.15;

% Constants
eSValue = 0.6; % eS
MxfITValue = 48; % MxfIT
KvVValue = 0.8*3e-5;% 8.7841e-06;% KvV KvG
KvGValue = 2*3e-5;
TauVValue = 1;
TauVAValue = 5;
TauBPValue = 1;
TauQValue = 10;
TauRValue = 1;
TauTAValue = 1;
TauITValue = 5;
TaupValue = 0.1;
RValue = 1.77e4; % R
VcValue = 19.2/2*1e-3;
VRValue = 133*1e-3;
VITValue = 0.05*1e-3;
s0Value = 5000/2;
kValue = 6000/2;
cpValue = 1000;
dAValue = 1.2;
MxDVAValue = 6.66;

% Derived
delta_hpITValue = ... % delta_hpIT
    (CoolProp.PropsSI('H','P',86.5e5,'S',...
    CoolProp.PropsSI('S','P',38e5,...
    'H',CoolProp.PropsSI('H','P',38e5,'Q',1,'CO2'),'CO2'),'CO2')-...
    CoolProp.PropsSI('H','P',38e5,'Q',1,'CO2'))/...
    (86.5e5-38e5);

save('constants.mat')