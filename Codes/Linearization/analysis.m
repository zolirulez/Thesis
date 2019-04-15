expA = sign(A).*exp(abs(A));
expA(A==0) = 0;
% ------------------------ MODAL ANALYSIS ---------------------------------
% Modal analysis
% Wt A = E Wt
% Wt A / Wt = E
% Dz = E z
% Dz = Wt A / Wt z
% Wt \ Dz = A / Wt z
% z = Wt x
[~,E,W] = eig(A);
Wt = W';
norm(Wt\E*Wt-A,Inf)
disp('The smaller this norm, the more reliable the modal matrix is')
% Separating real and imaginary values
Wt = real([Wt(1:6,:); real(Wt(7,:)); imag(Wt(7,:));...
    Wt(9,:); real(Wt(10,:)); imag(Wt(10,:)); Wt(12:end,:)]);
realexpeigA = sign(real(E)).*exp(abs(real(E)));
realexpeigA(real(E)==0) = 0;
imagexpeigA = sign(imag(E)).*exp(abs(imag(E)));
imagexpeigA(imag(E)==0) = 0;
% ----------------------- CONTROLLABILITY ---------------------------------
[AC,BC,CC,TC,KC] = ctrbf(A,B,C,1e-10);
% norm(AC-TC*A/TC,Inf)
% disp('The smaller this norm, the more reliable the controllability staircase form is')
rankCTRB = sum(KC)
disp(['If the system is not fully controllable, there are ' num2str(length(AC)-rankCTRB)...
    ' poles in the uncontrollable subspace'])
EnC = round(real(eig(AC(1:(length(AC)-rankCTRB),1:(length(AC)-rankCTRB)))),10);
disp(['There are ' num2str(sum(EnC==0))...
    ' integrators in the uncontrollable subspace, and ' num2str(sum(EnC<0))...
    ' a stable pole'])
expAC = sign(AC).*exp(abs(AC));
expAC(round(AC,10)==0) = 0;
GC = TC*G;
% ------------------------ OBSERVABILITY ----------------------------------
[AO,BO,CO,TO,KO] = obsvf(A,B,C,1e-10);
% norm(AO-TO*A/TO,Inf)
% disp('The smaller this norm, the more reliable the observability staircase form is')
rankOBSV = sum(KO)
disp(['If the system is not fully observable, there are ' num2str(length(AO)-rankOBSV)...
    ' poles in the unobservable subspace'])
EnO = round(real(eig(AO(1:(length(AO)-rankOBSV),1:(length(AO)-rankOBSV)))),10);
disp(['There are ' num2str(sum(EnO==0))...
    ' integrators in the unobservable subspace, and ' num2str(sum(EnO<0))...
    ' a stable pole'])
expAO = sign(AO).*exp(abs(AO));
expAO(round(AO,10)==0) = 0;
GO = TO*G;
% -------------------------- PLOTTING -------------------------------------
% System matrix
figure(1)
imagesc(expA,[-10 10])
colormap('jet')
colorbar
xlabel(char(x))
ylabel(['Diff' char(x)])
title('Pieceswise exponential of normalized system A, bounded by -10...10')
% Modal analysis
figure(2)
subplot(141)
imagesc(diag(imagexpeigA),[-10^0.5 10^0.5])
colormap('jet')
colorbar
xlabel('z')
ylabel('Dz')
title('imagexpeigA')
subplot(142)
imagesc(diag(realexpeigA),[-10 10])
colormap('jet')
colorbar
xlabel('z')
ylabel('Dz')
title('realexpeigA')
subplot(122)
imagesc(Wt,[-1 1])
colormap('jet')
colorbar
xlabel(char(x))
ylabel('z')
title('Wt, reordered')
% Controllability analysis
figure(3)
subplot(221)
imagesc(expAC,[-10^0.5 10^0.5])
colormap('jet')
colorbar
xlabel('z')
ylabel('Dz')
title('expAC')
subplot(223)
imagesc(CC,[-10^0 10^0])
colormap('jet')
colorbar
xlabel('z')
ylabel('y')
title('CC')
subplot(243)
imagesc(BC,[-10^-5 10^-5])
colormap('jet')
colorbar
xlabel(char(u))
ylabel('z')
title('BC')
subplot(244)
imagesc(GC,[-10^-5 10^-5])
colormap('jet')
colorbar
xlabel(char(d))
ylabel('z')
title('GdC')
% Observability analysis
figure(4)
subplot(221)
imagesc(expAO,[-10^0.5 10^0.5])
colormap('jet')
colorbar
xlabel('z')
ylabel('Dz')
title('expAO')
subplot(223)
imagesc(CO,[-10^1 10^1])
colormap('jet')
colorbar
xlabel('z')
ylabel('y')
title('CO')
subplot(243)
imagesc(BO,[-10^-5 10^-5])
colormap('jet')
colorbar
xlabel(char(u))
ylabel('z')
title('BO')
subplot(244)
imagesc(GO,[-10^-5 10^-5])
colormap('jet')
colorbar
xlabel(char(d))
ylabel('z')
title('GdO')