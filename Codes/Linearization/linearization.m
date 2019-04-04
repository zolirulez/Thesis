% ------------------------ LINEARIZATION ----------------------------------                              
% Jacobians
A = jacobian(Dx,x);
B = jacobian(Dx,u);
GBC = jacobian(Dx,BC);
GdDer = jacobian(Dx,dDer);
GdDis = jacobian(Dx,dDis);
Gd = jacobian(Dx,d);
C = jacobian(y,x);
D = jacobian(y,u);
EBC = jacobian(y,BC);
EdDer = jacobian(y,dDer);
EdDis = jacobian(y,dDis);
Ed = jacobian(y,d);
% ------------------------ SUBSTITUTIONS ----------------------------------
mx4sub = struct('A',A,'B',B,'C',C,'Gd',Gd);
fields = fieldnames(mx4sub);
substitution;
A = double(mx4sub.A);
B = double(mx4sub.B);
Gd = double(mx4sub.Gd);
C = double(mx4sub.C);
D = double(D);
% ------------------------- NORMALIZING -----------------------------------
% Normalizing with maximum required deviations
% x = [DVA; p1; h1; d1; T1; TA2; Dm21; p2; h2; d2; T2; TA1; BP; DmV; pR; hR; dR; DmG; DmL; delta_hJ; delta_h2];
condA = cond(A)
disp('The condition of the matrix indicates the relative sensitivies within the system')
DVBound = 1;
pBound = 5*10^5;
hBound = 20*10^3;
BPBound = 0.1;
dBound = 10;
TBound = 5;
DmBound = 0.1;
T = diag(1./[DVBound;...
    pBound; hBound; dBound; TBound; TBound;...
    DmBound; ...
    pBound; hBound; dBound; TBound; TBound;...
    BPBound; DmBound; ...
    pBound; hBound; dBound;...
    DmBound; DmBound;...
    hBound; hBound]);
% z = Tx
% Dz = TA/Tx
A = T*A/T;
B = T*B;
GBC = T*GBC;
GdDer = T*GdDer;
GdDis = T*GdDis;
Gd = T*Gd;
C = C/T;
condA = cond(A)
disp('The condition of the matrix decreased by normalizing')
