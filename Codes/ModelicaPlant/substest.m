clearvars
% Part 1: can be done in advance
addpath('C:\Users\u375749\Documents\Thesis\Codes\Linearization')
addpath('C:\Users\u375749\Documents\Thesis\Codes\MatlabPlant')
syms A B C D
AValue = 1;
BValue = 2;
CValue = 4;
DValue = 5;
% modelcreation_simplified
ABCD = [A B; C D];
symvars = symvar(ABCD);
symvarsstr = ['[', strjoin(arrayfun(@char, symvars, 'uniform', 0),', '), ']'];
endindex = find(any([symvarsstr == ','; symvarsstr == ']']));
startindex = find(any([symvarsstr == '['; symvarsstr == ' ']));
value = cell(length(symvars),1);
for it = 1:length(symvars)
    value{it} = [symvarsstr(startindex(it)+1:endindex(it)-1) 'Value'];
end
% Part 2: done online
subsarray = [];
for it = 1:length(symvars)
    subsarray = [subsarray eval(value{it})];
end
double(subs(ABCD,symvars,subsarray))


