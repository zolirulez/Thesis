function ABCD = LTVsystemDescription(u,x)

constants = load('constants.mat')
fields = fieldnames(constants);
for i=1:length(fields)
    eval([fields{i} '=constants.' fields{i} ]);
end

system = struct('A',A,'B',B,'C',C,'D',D);
fields = fieldnames(system);

% Input values
BPRValue = u(1);
CRVValue = u(2);
CRAValue = u(3);
CRITValue = u(4);
delta_hHRValue = u(5);
dBPValue = u(6);
dGValue = u(7);
hGValue = u(8);
hLValue = u(9);
TA0Value = u(10);
hMTValue = u(11);

% State values
DVAValue = x(1);
p1Value = x(2);
h1Value = x(3);
d1Value = x(4);
TA2Value = x(5);
Dm21Value = x(6);
p2Value = x(7);
h2Value = x(8);
d2Value = x(9);
TA1Value = x(10);
BPValue = x(11);
DmVValue = x(12);
pRValue = x(13);
hRValue = x(14);
dRValue = x(15);
DmGValue = x(16);
delta_h2Value = x(17);
DmQValue = x(18);

% Derived values
delta_hpITValue = ... % delta_hpIT
    (CoolProp.PropsSI('H','P',p1Value,'S',...
    CoolProp.PropsSI('S','P',pRValue,'H',hGValue,'CO2'),'CO2')-hGValue)/...
    (p1Value-pRValue);

% Substitutions
subsarray = [];
for it = 1:length(symvars)
    subsarray = [subsarray eval(value{it})];
end
for it = 1:numel(fields)
    system.(fields{it}) = double(subs(system.(fields{it}),symvars,subsarray));
end

nx = length(A);
nu = size(B,2);
ny = size(C,1);
ABCD = [reshape(system.A,nx*nx,1); ...
    reshape(system.B,nx*nu,1); ...
    reshape(system.C,ny*nx,1); ...
    reshape(system.D,ny*nu,1)];