figure, plot(start+2001:length(outcorrecord)+start+2000,outcorrecord')
hold on
plot(start:finish,outrecord')
plot([detectiontime-300 detectiontime-300],[0 1.5e5],'--')
hold off