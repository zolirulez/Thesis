function xp = constrainer(x,ex,k)
% Steady state constrainer function.
% Input x, limit ex, double of the maximum deviation: k
% Zoltan Mark Pinter, Master Thesis, 2019

if x > ex
    beta = 1/k*(x-ex);
    xp = (1-beta)*x + beta*ex;
    xp = max(xp,ex);
    warning('Constraint applied')
else
    xp = x;
end