% RLS
rlsInitial.t = ones(6,1)*[1 1e-3];  % 10
rlsInitial.P = diag(max(rlsInitial.t')')/10; 
rls = RecursiveLeastSquares;
lambda = 1 - 1e-4;
weight = eye(size(rlsInitial.t,2));
rls.initialize(lambda,weight,rlsInitial);
nr = size(rlsInitial.t,1);
npy = size(rlsInitial.t,2);

% % RLS2
% rlsInitial2.t = [5000; 6000; 1e3; 1e3; 1e4; 1e4; 1e3; 600; 1e2; 1e2; 1e3; 1e3; 1e2]*[1 1e-3];  % 10
% rlsInitial2.P = diag(max(rlsInitial2.t')')/10; 
% rls2 = RecursiveLeastSquares;
% lambda2 = 1;
% weight2 = eye(size(rlsInitial2.t,2));
% rls2.initialize(lambda2,weight2,rlsInitial2);