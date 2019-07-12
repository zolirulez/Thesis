function xp = constrainer(x,ex,k)

if x > ex
    beta = 1/k*(x-ex);
    xp = (1-beta)*x + beta*ex;
    xp = max(xp,ex);
else
    xp = x;
end