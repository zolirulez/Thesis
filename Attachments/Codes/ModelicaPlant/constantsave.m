function [A,B,C,D] = constantsave(A,B,C,D,Q)
% Function to speed up later function substitutions
% Zoltan Mark Pinter, Master Thesis, 2019

load transferparam
% Cholesky factorisation
Lq = chol(Q,'lower');

% Dimensions
nx = length(A);
nu = size(B,2);
ny = size(C,1);

ABCD = [A B; C D];

% Constants
c.eSValue = 0.65;
c.MxfITValue = 48; 
c.KvVValue = 0.8*3e-5; 
c.KvGValue = 2*3e-5;
c.TauAValue = 150;
c.TaupValue = 1;
c.TauhValue = 5;
c.VcValue = 19.2*1e-3*2;
c.VRValue = 133*1e-3*2;
c.V_ITValue = 0.05*1e-3*2; % Multiplication by 2: described in substitution-simplified.m
c.s0Value = 5000;
c.kValue = 6000;
c.cpValue = 1000;
c.dAValue = 1.25;
c.MxDVAValue = 5.2;

% Substitutions
fields = fieldnames(c);
for it = 1:numel(fields)
    eval(['syms ' fields{it}(1:end-5)]);
    ABCD = subs(ABCD,eval(['{' num2str(fields{it}(1:end-5)) '}']),c.(fields{it}));
end

% Converting symbolic function to Matlab function
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

% Save constants for speeding up LTVsystemDescription.m
save('constants.mat')