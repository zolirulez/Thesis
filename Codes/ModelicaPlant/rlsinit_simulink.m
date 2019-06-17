% RLS
rlsInitial.t = [5000; 6000; 1e3; 1e3; 1e3]*[1 1e-3];  % 10
rlsInitial.P = diag(max(rlsInitial.t')')/10; 
rls = RecursiveLeastSquares;
lambda = 1 - 1e-4;
weight = eye(size(rlsInitial.t,2));
rls.initialize(lambda,weight,rlsInitial);