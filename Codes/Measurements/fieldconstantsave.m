function [A,B,C,D] = fieldconstantsave(A,B,C,D,Q)

load transferparam

% c.delta_ph1Value = 13;
% c.delta_ph2Value = 17;
% c.delta_phRValue = 8;
% c.delta_pd1Value = 1.25*1e4;
% c.delta_pd2Value = 0.46*1e4;
% c.delta_pdRValue = 0.50*1e4;
% c.delta_Th1Value = 6.6*1e-4;
% c.delta_Th2Value = 7.0*1e-4;
% c.delta_Td1Value = 0.16;
% c.delta_Td2Value = 0.15;

Lq = chol(Q,'lower');

% Dimensions
nx = length(A);
nu = size(B,2);
ny = size(C,1);

ABCD = [A B; C D];

% Constants
c.eSValue = 0.65; % eS
c.MxfITValue = 48; % MxfIT
c.KvVValue = 0.8*3e-5; 
c.KvGValue = 1.6*3e-5; 
c.TauTAValue = 1;
c.TaupValue = 0.1;
c.TauhValue = 5;
c.VcValue = 19.2*1e-3;
c.VRValue = 115*1e-3;
c.V_ITValue = 3.3/3600/48*2; % Multiplication by 2: described in substitution-simplified.m
c.s0Value = 5000;
c.kValue = 6000;
c.cpValue = 1000;
c.dAValue = 1.25;
c.MxDVAValue = 5.2;

% Derived
% c.delta_hpITValue = ... % delta_hpIT
%     (CoolProp.PropsSI('H','P',86.5e5,'S',...
%     CoolProp.PropsSI('S','P',38e5,...
%     'H',CoolProp.PropsSI('H','P',38e5,'Q',1,'CO2'),'CO2'),'CO2')-...
%     CoolProp.PropsSI('H','P',38e5,'Q',1,'CO2'))/...
%     (86.5e5-38e5);

fields = fieldnames(c);
for it = 1:numel(fields)
    eval(['syms ' fields{it}(1:end-5)]);
    ABCD = subs(ABCD,eval(['{' num2str(fields{it}(1:end-5)) '}']),c.(fields{it}));
end

symvars = symvar(ABCD);
symvarsstr = ['(', strjoin(arrayfun(@char, symvars, 'uniform', 0),', '), ')'];
endindex = find(any([symvarsstr == ','; symvarsstr == ')']));
startindex = find(any([symvarsstr == '('; symvarsstr == ' ']));
value = [];
for it = 1:length(symvars)
    value = [value symvarsstr(startindex(it)+1:endindex(it)-1) 'Value,'];
end
value = value(1:end-1);

ABCD = matlabFunction(ABCD);

save('constants.mat')