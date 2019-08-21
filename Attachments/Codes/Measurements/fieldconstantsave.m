function [A,B,C,D] = fieldconstantsave(A,B,C,D,Q)
% Function to speed up later function substitutions

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
c.KvGValue = 1.6*3e-5; 
c.TauAValue = 1;
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