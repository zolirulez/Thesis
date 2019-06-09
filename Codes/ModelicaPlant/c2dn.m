function [Ad,Bd,Qd]=c2dn(Ac,Bc,Gc,Ts)
[nx,nu]=size(Bc);
M = [Ac Bc; zeros(nu,nx+nu)];
if ~any(isnumeric(M))
    disp('Untracked variable(s):')
    symvar(M)
    error('Symbolic matrix for exponential matrix calculation')
end
Phi = expm(M*Ts);
Ad = Phi(1:nx,1:nx);
Bd = Phi(1:nx,nx+1:nx+nu);

% M = [-Ac' Gc*Gc'; zeros(nx,nx) Ac];
% Phi = expm(M*Ts);
% Abar = Phi(nx+1:nx+nx,nx+1:nx+nx);
% Qd = Ad'*Phi(1:nx,nx+1:nx+nx);
Qd = Gc*Gc';