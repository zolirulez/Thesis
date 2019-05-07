t = start:Ts:finish;
figure(1)
clf
subplot(321)
hold on
plot(t,record(2,:)/1e5)
plot(t,record(7,:)/1e5)
hold off
ylim([30 100])
xlabel('Time [s]')
ylabel('Pressure [bar]')
subplot(322)
hold on
plot(t,record(3,:)/1000)
plot(t,record(8,:)/1000)
plot(t,record(11,:)/1000)
plot(t,record(13,:)/1000)
hold off
ylim([0 600])
xlabel('Time [s]')
ylabel('Enthalpy [kJ/kg]')
subplot(323)
hold on
plot(t,record(5,:)-273.15)
hold off
ylim([-10 120])
xlabel('Time [s]')
ylabel('Temperature [C]')
subplot(324)
hold on
plot(t,record(6,:))
plot(t,record(10,:))
plot(t,record(12,:))
hold off
ylim([0 0.5])
xlabel('Time [s]')
ylabel('Mass flow rate [kg/s]')
subplot(325)
hold on
plot(t,record(4,:))
plot(t,record(9,:))
hold off
ylim([0 900])
xlabel('Time [s]')
ylabel('Density [kg/m^3]')
subplot(326)
hold on
plot(t,record(1,:))
hold off

figure(2)
subplot(311)
plot(t,record(nx+1:nx+nx,:))
xlabel('Time [s]')
ylabel('Eigenvalues of P_1')
subplot(312)
plot(t,record(nx+nx+1:nx+nx+nx,:)./diag(noise.Q))
xlabel('Time [s]')
ylabel('Relative state correction of K_xe')
legend('DVA','p1','h1','d1','TA1','DmV','pR','hR','dR','DmG','delta_h2','DmQ','hHR');
subplot(3,5,11)
plot(t,record(nx+nx+nx+1,:)/1e5,t,record(nx+nx+nx+3,:)/1e5)
legend('p_2','p_R')
xlabel('Time [s]')
ylabel('Innovation [bar]')
subplot(3,5,12)
plot(t,record(nx+nx+nx+2,:)/1e3,t,record(nx+nx+nx+4,:)/1e3,t,record(nx+nx+nx+5,:)/1e3)
legend('h_B_P''h_R','h_H_R')
xlabel('Time [s]')
ylabel('Innovation [kJ/kg]')
subplot(3,5,13)
plot(t,record(nx+nx+nx+6,:))
legend('T_A_1')
xlabel('Time [s]')
ylabel('Innovation [C]')
subplot(3,5,14)
plot(t,record(nx+nx+nx+7,:))
legend('Dm_Q')
xlabel('Time [s]')
ylabel('Innovation [kg/s]')
subplot(3,5,15)
plot(t,record(nx+nx+nx+8,:))
legend('d_R')
xlabel('Time [s]')
ylabel('Innovation [kg/m^3]')