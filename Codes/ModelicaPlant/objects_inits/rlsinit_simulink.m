% RLS
rlsInitial.t = ones(7,1)*[1 1e-3];  % 10
rlsInitial.P = diag(max(rlsInitial.t')')/10; 
rls = RecursiveLeastSquares;
lambda = 1 - 1e-4;
weight = eye(size(rlsInitial.t,2));
rls.initialize(lambda,weight,rlsInitial);
nr = size(rlsInitial.t,1);
npy = size(rlsInitial.t,2);