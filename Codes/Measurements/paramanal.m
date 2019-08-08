function figurehandle = paramanal(beststr,Z,direction)

ns8=beststr(8-2,:) %[3 4 1 0 0 1]
th8=pem(Z,ns8);
res8=pe(Z,th8);
[~,p8] = estpres(th8); % might be referred
cond8 = cond(p8);
ns7=beststr(7-2,:) %[3 2 1 0 0 1]
th7=pem(Z,ns7);
res7=pe(Z,th7);
[~,p7] = estpres(th7); % might be referred
cond7 = cond(p7);
ns6=beststr(6-2,:) %[3 2 1 0 0 1]
th6=pem(Z,ns6);
res6=pe(Z,th6);
[~,p6] = estpres(th6) % might be referred
cond6 = cond(p6);
ns5=beststr(5-2,:) %[2 2 1 0 0 1]
th5=pem(Z,ns5);
res5=pe(Z,th5);
[~,p5] = estpres(th5) % might be referred
cond5 = cond(p5);
cond85 = [cond8; cond7; cond6; cond5]

h=figure(6)
subplot(221)
resid(Z,th8);
title('Residue corr., 8 param.')
subplot(222)
resid(Z,th7);
title('Residue corr., 7 param.')
subplot(223)
resid(Z,th6);
title('Residue corr., 6 param.')
subplot(224)
resid(Z,th5);
title('Residue corr., 5 param.')
if direction == 2
    % saveas(h,'autocorr.png')
end
figurehandle = figure(30);
% h7=figure(7)
% subplot(221)
% hold on
% zpplot(th2zp(th8,0),2.58);
% h = findobj(gca,'type','line');
% if direction == 2
%     set(h([floor(end/2):end-2]),'Color','r');
%     title('From noise to output, 8 param.')
% end
% hold off
% subplot(222)
% hold on
% zpplot(th2zp(th7,0),2.58);
% h = findobj(gca,'type','line');
% if direction == 2
%     set(h([floor(end/2):end-2]),'Color','r');
%     title('From noise to output, 7 param.')
% end
% hold off
subplot(235)
hold on
zpplot(th2zp(th6,0),2.58);
h = findobj(gca,'type','line');
if direction == 2
    set(h([floor(end/2):end-2]),'Color','r');
    title('From noise to output, 6 param.')
end
hold off
xlabel('Real [rad^-^1]')
ylabel('Imaginary [i rad^-^1]')
subplot(236)
hold on
zpplot(th2zp(th5,0),2.58);
h = findobj(gca,'type','line');
if direction == 2
    set(h([floor(end/2):end-2]),'Color','r');
    title('From noise to output, 5 param.')
    % saveas(h7,'noisetoouput.png')
end
hold off
xlabel('Real [rad^-^1]')
ylabel('Imaginary [i rad^-^1]')
% h8=figure(8)
% subplot(221)
% hold on
% zpplot(th2zp(th8,1),2.58);
% h = findobj(gca,'type','line');
% if direction == 2
%     set(h([floor(end/2):end-2]),'Color','r');
%     title('From input to output, 8 param.')
% end
% hold off
% subplot(222)
% hold on
% zpplot(th2zp(th7,1),2.58);
% h = findobj(gca,'type','line');
% if direction == 2
%     set(h([floor(end/2):end-2]),'Color','r');
%     title('From input to output, 7 param.')
% end
% hold off
subplot(232)
hold on
zpplot(th2zp(th6,1),2.58);
h = findobj(gca,'type','line');
if direction == 2
    set(h([floor(end/2):end-2]),'Color','r');
    title('From input to output, 6 param.')
end
hold off
xlabel('Real [rad^-^1]')
ylabel('Imaginary [i rad^-^1]')
subplot(233)
hold on
zpplot(th2zp(th5,1),2.58);
h = findobj(gca,'type','line');
if direction == 2
    set(h([floor(end/2):end-2]),'Color','r');
    title('From input to output, 5 param.')
    % saveas(h8,'inputtoouput.png')
end
hold off
xlabel('Real [rad^-^1]')
ylabel('Imaginary [i rad^-^1]')
h=figure(40)
hold on
subplot(221)
hold on
[pxx,f,pxxc] = periodogram(res8, 'ConfidenceLevel', 0.99);
cpxx = cumsum(pxx)./sum(pxx);
plot(f,cpxx)
plot(f,f./max(f)+0.05,'r--')
plot(f+0.05*max(f),f./max(f),'r--')
hold off
title('Cum. Per. of residual 8')
xlabel('Frequency')
ylabel('Power')
legend('Forward test','Backward test')
subplot(222)
hold on
[pxx,f,pxxc] = periodogram(res7, 'ConfidenceLevel', 0.99);
cpxx = cumsum(pxx)./sum(pxx);
plot(f,cpxx)
plot(f,f./max(f)+0.05,'r--')
plot(f+0.05*max(f),f./max(f),'r--')
hold off
title('Cum. Per. of residual 7')
xlabel('Frequency')
ylabel('Power')
legend('Forward test','Backward test')
subplot(223)
hold on
[pxx,f,pxxc] = periodogram(res6, 'ConfidenceLevel', 0.99);
cpxx = cumsum(pxx)./sum(pxx);
plot(f,cpxx)
plot(f,f./max(f)+0.05,'r--')
plot(f+0.05*max(f),f./max(f),'r--')
hold off
title('Cum. Per. of residual 6')
xlabel('Frequency')
ylabel('Power')
legend('Forward test','Backward test')
subplot(224)
hold on
[pxx,f,pxxc] = periodogram(res5, 'ConfidenceLevel', 0.99);
cpxx = cumsum(pxx)./sum(pxx);
plot(f,cpxx)
plot(f,f./max(f)+0.05,'r--')
plot(f+0.05*max(f),f./max(f),'r--')
hold off
title('Cum. Per. of residual 5')
xlabel('Frequency')
ylabel('Power')
legend('Forward test','Backward test')
if direction == 2
    % saveas(h,'periodogram.png')
end

% figure(40)
% hold on
% subplot(221)
% probplot(res8);
% subplot(222)
% probplot(res7);
% subplot(223)
% probplot(res6);
% subplot(224)
% probplot(res5);
% hold off