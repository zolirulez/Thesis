t = start:Ts:finish;
figure(1)
clf
subplot(321)
hold on
plot(t,record(1,:)/1e5)
plot(t,record(5,:)/1e5)
hold off
ylim([30 100])
xlabel('Time [s]')
ylabel('Pressure [bar]')
subplot(322)
hold on
plot(t,record(2,:)/1000)
plot(t,record(6,:)/1000)
plot(t,record(8,:)/1000)
hold off
ylim([0 600])
xlabel('Time [s]')
ylabel('Enthalpy [kJ/kg]')
subplot(323)
hold on
plot(t,record(4,:)-273.15)
hold off
ylim([-10 120])
xlabel('Time [s]')
ylabel('Temperature [C]')
subplot(324)
hold on
plot(t,record(9,:))
hold off
ylim([0 0.5])
xlabel('Time [s]')
ylabel('Mass flow rate [kg/s]')
subplot(325)
hold on
plot(t,record(3,:))
plot(t,record(7,:))
hold off
ylim([0 900])
xlabel('Time [s]')
ylabel('Density [kg/m^3]')

figure(2)
subplot(311)
plot(t,record(nx+1:nx+nx,:))
xlabel('Time [s]')
ylabel('Eigenvalues of P_1')
subplot(312)
plot(t,record(nx+nx+1:nx+nx+nx,:))
xlabel('Time [s]')
ylabel('State correction of K_xe')
legend('p1','h1','d1','TA1','pR','hR','dR','delta_h','DmQ');
subplot(3,4,9)
plot(t,record(nx+nx+nx+1,:)/1e5,t,record(nx+nx+nx+3,:)/1e5)
legend('p_1','p_R')
xlabel('Time [s]')
ylabel('Innovation [bar]')
subplot(3,4,10)
plot(t,record(nx+nx+nx+2,:)/1e3,t,record(nx+nx+nx+4,:)/1e3)
legend('h_B_P','h_R')
xlabel('Time [s]')
ylabel('Innovation [kJ/kg]')
subplot(3,4,11)
plot(t,record(nx+nx+nx+5,:))
legend('Dm_Q')
xlabel('Time [s]')
ylabel('Innovation [kg/s]')
subplot(3,4,12)
plot(t,record(nx+nx+nx+6,:))
legend('d_R')
xlabel('Time [s]')
ylabel('Innovation [kg/m^3]')

figure(3)
plot(t,record(nx+nx+nx+ny+1:nx+nx+nx+ny+2,:))
legend('s_0','k')
xlabel('Time [s]')
ylabel('Parameters')