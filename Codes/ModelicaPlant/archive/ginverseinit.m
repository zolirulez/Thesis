function ginv = ginverseinit(g)

load constants

syms ginv
for it = 1:length(g)
    ginv(it,1) = finverse(g(it));
end
fields = fieldnames(c);
for it = 1:numel(fields)
    eval(['syms ' fields{it}(1:end-5)]);
    ginv = subs(ginv,eval(['{' num2str(fields{it}(1:end-5)) '}']),c.(fields{it}));
end

