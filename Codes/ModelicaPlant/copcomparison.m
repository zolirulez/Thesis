% figure,
% load faultignore_long
% plot(COP_sim.signals.values(:,1)')
% hold on
% load faultcontrol_nobatch_long
% plot(COP_sim.signals.values(:,1)')
% load faultcontrol_long
% plot(COP_sim.signals.values(:,1)')
% hold off

figure,
load faultignore_long
plot(residual_sim.signals.values(:,2)')
hold on
load faultcontrol_nobatch_long
plot(residual_sim.signals.values(:,2)')
load faultcontrol_long
plot(residual_sim.signals.values(:,2)')
hold off