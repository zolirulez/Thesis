x0 = 0:1500;
h = figure(1)
set(h, 'Position',  [100, 100, 100+300, 100+200])
plot(x0,1000*ones(length(x0),1),'--')
x1 = 0:1000;
x2 = 1000:1225;
x3 = 1225:1500;
hold on
plot(x1,x1,'r')
a = 1/50000*(x2-1000).^2;
plot(x2,(1-a).*x2+a*1000,'r')
plot(x3,1000*ones(length(x3),1),'r')
hold off
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
xlabel('$X_{raw}$','Interpreter','latex')
ylabel('$X_{f}$','Interpreter','latex')
saveas(h,'constrainer.png')

h = figure(2);
set(h, 'Position',  [100, 100, 100+600, 100+200])
x = [-0.5:.01:0.5];
y = normpdf(x,0,(1/3)^2);
plot(x,y)
y = normpdf(x,0,(1/3)^2*2);
hold on
plot(x,y,'r--')
y = normpdf(x,0.2,(1/3)^2*2);
plot(x,y,'r-')
plot([0 0],[0 5],'g--')
plot([0.2 0.2],[0 5],'g-')
hold off
ylim([0 4])
xticks([])
yticks([])
legend({'Normal operation cluster','Faulty operation cluster, no fault',...
    'Faulty operation cluster, fault','Normal operation residual mean'...
    ,'Faulty operation residual mean'},'Location','northwest')
xlabel('Residual value')
ylabel('Probability density')
saveas(h,'clusters.png')