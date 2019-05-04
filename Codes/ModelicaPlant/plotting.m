t = start:Ts:finish;
figure(1)
clf
subplot(321)
hold on
plot(t,record(2,:)/1e5)
plot(t,record(7,:)/1e5)
plot(t,record(13,:)/1e5)
hold off
ylim([30 100])
xlabel('Time [s]')
ylabel('Pressure [bar]')
subplot(322)
hold on
plot(t,record(3,:)/1000)
plot(t,record(8,:)/1000)
plot(t,record(14,:)/1000)
plot(t,record(17,:)/1000)
hold off
ylim([200 600])
xlabel('Time [s]')
ylabel('Enthalpy [kJ/kg]')
subplot(323)
hold on
plot(t,record(5,:)-273.15)
plot(t,record(10,:)-273.15)
hold off
ylim([-10 120])
xlabel('Time [s]')
ylabel('Temperature [C]')
subplot(324)
hold on
plot(t,record(6,:))
plot(t,record(12,:))
plot(t,record(16,:))
plot(t,record(18,:))
hold off
ylim([0 0.5])
xlabel('Time [s]')
ylabel('Mass flow rate [kg/s]')
subplot(325)
hold on
plot(t,record(4,:))
plot(t,record(9,:))
plot(t,record(15,:))
hold off
ylim([0 900])
xlabel('Time [s]')
ylabel('Density [kg/m^3]')
subplot(326)
hold on
plot(t,record(1,:))
plot(t,record(11,:))
hold off

figure(2)
subplot(311)
plot(t,record(nx+1:nx+nx,:))
xlabel('Time [s]')
ylabel('Eigenvalues of P_1')
subplot(312)
plot(t,record(nx+nx+1:nx+nx+nx,:))
ylim([-3e3 3e3])
xlabel('Time [s]')
ylabel('State correction of K_xe')
legend('DVA','p1','h1','d1','TA2','Dm21','p2','h2','d2','TA1','BP','DmV','pR','hR','dR','DmG','delta_h2','DmQ');
subplot(313)
plot(t,record(nx+nx+nx+1:nx+nx+nx+ny,:))
ylim([-3e3 3e3])
legend('p_2','h_B_P','p_R','h_R','h_H_R','T_A_1','Dm_Q','Dm_2_1')
xlabel('Time [s]')
ylabel('Innovation')