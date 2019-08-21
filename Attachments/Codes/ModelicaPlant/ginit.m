function g = ginit(g)
% Function to initialize nonlinear output function g

load constants
fields = fieldnames(c);
for it = 1:numel(fields)
    eval(['syms ' fields{it}(1:end-5)]);
    g = subs(g,eval(['{' num2str(fields{it}(1:end-5)) '}']),c.(fields{it}));
end

