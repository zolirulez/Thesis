% Script to make LTI analysis (stability,ctrb,obsv) at the picked set point
% Zoltan Mark Pinter, Master Thesis, 2019


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
Wt = real([Wt(1:4,:); real(Wt(5,:)); imag(Wt(5,:)); Wt(7:end,:)]);
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
%GC = TC*G;
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
%GO = TO*G;
% -------------------------- PLOTTING -------------------------------------
% System matrix
handle = figure(1);
set(handle, 'Position',  [100, 100, 100+500, 100+300])
imagesc(expA,[-10 10])
colormap('jet')
colorbar
xticks(1:8)
set(gca,'TickLabelInterpreter','latex')
xticklabels({'$p_{GC}$','$h_{GC}$','$\rho_{GC}$',...
    '$T_{A,GC}$','$p_{R}$','$h_{R}$','$\rho_{GC}$',....
    '$\Delta h$'})
yticks(1:8)
set(gca,'TickLabelInterpreter','latex')
yticklabels({'$\dot{p}_{GC}$','$\dot{h}_{GC}$','$\dot{\rho}_{GC}$',...
    '$\dot{T}_{A,GC}$','$\dot{p}_{R}$','$\dot{p}_{R}$','$\dot{\rho}_{GC}$',....
    '$\dot{\Delta h}$'})
title('Piecewise exponential of normalized system matrix')
saveas(handle,'amatrix.png')
% Modal analysis
handle = figure(2);
set(handle, 'Position',  [100, 100, 100+1200, 100+300])
subplot(141)
imagesc(diag(imagexpeigA),[-10^0.5 10^0.5])
colormap('jet')
colorbar
yticks([])
xlabel('\textbf{z}','Interpreter','latex')
ylabel('$\mathbf{\dot{z}}$','Interpreter','latex')
title('Imaginary values')
subplot(142)
imagesc(diag(realexpeigA),[-10 10])
colormap('jet')
colorbar
yticks([])
xlabel('\textbf{z}','Interpreter','latex')
ylabel('$\mathbf{\dot{z}}$','Interpreter','latex')
title('Real values')
subplot(122)
imagesc(Wt,[-1 1])
xticks(1:8)
set(gca,'TickLabelInterpreter','latex')
xticklabels({'$p_{GC}$','$h_{GC}$','$\rho_{GC}$',...
    '$T_{A,GC}$','$p_{R}$','$h_{R}$','$\rho_{GC}$',....
    '$\Delta h$'})
yticks([])
colormap('jet')
colorbar
ylabel('$\mathbf{z}$','Interpreter','latex')
title('W^T')
saveas(handle,'eigen.png')
% Controllability analysis
handle = figure(3);
set(handle, 'Position',  [100, 100, 100+800, 100+300])
subplot(221)
imagesc(expAC,[-10^0.5 10^0.5])
colormap('jet')
colorbar
xlabel('$\mathbf{z}$','Interpreter','latex')
ylabel('$\mathbf{\dot{z}}$','Interpreter','latex')
title('System matrix A_C')
subplot(223)
imagesc(CC,[-10^0 10^0])
colormap('jet')
colorbar
xticks([])
yticks([])
xlabel('$\mathbf{z}$','Interpreter','latex')
ylabel('$\mathbf{y}$','Interpreter','latex')
title('C_C')
subplot(243)
imagesc(BC,[-10^-5 10^-5])
colormap('jet')
colorbar
xticks([])
yticks([])
xlabel('$\mathbf{u}$','Interpreter','latex')
ylabel('$\mathbf{z}$','Interpreter','latex')
title('B_C')
saveas(handle,'ctrb.png')
% Observability analysis
handle = figure(4);
set(handle, 'Position',  [100, 100, 100+800, 100+300])
subplot(221)
imagesc(expAO,[-10^0.5 10^0.5])
colormap('jet')
colorbar
xticks([])
yticks([])
xlabel('$\mathbf{z}$','Interpreter','latex')
ylabel('$\mathbf{\dot{z}}$','Interpreter','latex')
title('System matrix A_O')
subplot(223)
imagesc(CO,[-10^1 10^1])
colormap('jet')
colorbar
xticks([])
yticks([])
xlabel('$\mathbf{z}$','Interpreter','latex')
ylabel('$\mathbf{y}$','Interpreter','latex')
title('C_O')
subplot(243)
imagesc(BO,[-10^-5 10^-5])
colormap('jet')
colorbar
xticks([])
yticks([])
xlabel('$\mathbf{u}$','Interpreter','latex')
ylabel('$\mathbf{z}$','Interpreter','latex')
title('B_O')
saveas(handle,'obsv.png')