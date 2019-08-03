figure,
load faultignore_long2_ta0
plot(COP_sim.signals.values(:,1)')
hold on
load faultcontrol_nobatch_long2_ta0
plot(COP_sim.signals.values(:,1)')
load faultcontrol_long2_ta0
plot(COP_sim.signals.values(:,1)')
hold off

figure,
load faultignore_long2_ta0
subplot(311)
plot(squeeze(residual_sim.signals.values)')
load faultcontrol_nobatch_long2_ta0
(var(squeeze(residual_sim.signals.values(:,:,1000:5000))'))
subplot(312)
plot(squeeze(residual_sim.signals.values)')
load faultcontrol_long2_ta0
(var(squeeze(residual_sim.signals.values(:,:,1000:5000))'))
subplot(313)
plot(squeeze(residual_sim.signals.values)')