% ------------------------- NORMALIZING -----------------------------------
% Normalizing with maximum required deviations
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