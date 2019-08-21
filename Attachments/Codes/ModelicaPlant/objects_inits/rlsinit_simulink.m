% RLS object initialization

% Parameter vector
rlsInitial.t = ones(7,1)*[1 1e-3];  
% Covariance
rlsInitial.P = diag(max(rlsInitial.t')')/10; 
% Object instance
rls = RecursiveLeastSquares;
% Forgetting factor
lambda = 1 - 1e-4;
% Initialization
rls.initialize(lambda,rlsInitial);
% Dimensions
nr = size(rlsInitial.t,1);
npy = size(rlsInitial.t,2);