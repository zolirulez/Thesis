% ------------------------- NORMALIZING -----------------------------------
% Normalizing with maximum required deviations
% x = [DVA; p1; h1; d1; TA1; BP; DmV; pR; hR; dR; DmG; delta_h; DmQ]
condA = cond(A)
disp('The condition of the matrix indicates the relative sensitivies within the system')
T = diag(1./[DVBound;...
    pBound; hBound; dBound; TBound;...
    BPBound; DmBound; ...
    pBound; hBound; dBound;...
    DmBound;...
    hBound;  DmBound]);
% z = Tx
% Dz = TA/Tx
A = T*A/T;
B = T*B;
G = T*G;
C = C/T;
condA = cond(A)
disp('The condition of the matrix decreased by normalizing')