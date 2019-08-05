figure,
load faultignore
plot(COP_sim.signals.values(:,1)')
hold on
load faultcontrol_nobatch
plot(COP_sim.signals.values(:,1)')
load faultcontrol
plot(COP_sim.signals.values(:,1)')
hold off

figure,
load faultignore
subplot(311)
plot(squeeze(residual_sim.signals.values)')
load faultcontrol_nobatch
(var(squeeze(residual_sim.signals.values(:,:,1000:5000))'))
subplot(312)
plot(squeeze(residual_sim.signals.values)')
load faultcontrol
(var(squeeze(residual_sim.signals.values(:,:,1000:5000))'))
subplot(313)
plot(squeeze(residual_sim.signals.values)')