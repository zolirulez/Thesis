% ------------------------- NORMALIZING -----------------------------------
% Normalizing with maximum required deviations
condA = cond(A)
disp('The condition of the matrix indicates the relative sensitivies within the system')
T = diag(1./diag(Qcont));
% z = Tx
% Dz = TA/Tx
A = T*A/T;
B = T*B;
C = C/T;
condA = cond(A)
disp('The condition of the matrix decreased by normalizing')