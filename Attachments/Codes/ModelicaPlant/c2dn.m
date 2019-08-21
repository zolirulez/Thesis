function [Ad,Bd,Qd]=c2dn(Ac,Bc,Gc,Ts)
% Function to convert continuous time state space model to discrete time,
% possibly with covariance
% Prepared for the Model Predictive Control course at DTU

% System matrices
[nx,nu]=size(Bc);
M = [Ac Bc; zeros(nu,nx+nu)];
Phi = expm(M*Ts);
Ad = Phi(1:nx,1:nx);
Bd = Phi(1:nx,nx+1:nx+nu);

% ------- Covariance --------

% First case: the text book way
% M = [-Ac' Gc*Gc'; zeros(nx,nx) Ac];
% Phi = expm(M*Ts);
% Abar = Phi(nx+1:nx+nx,nx+1:nx+nx);
% Qd = Ad'*Phi(1:nx,nx+1:nx+nx);

% Second case: simplification
Qd = Gc*Gc';