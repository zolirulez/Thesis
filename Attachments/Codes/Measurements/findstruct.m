function [AIC,BIC,FPE,beststrt,figurehandle] = findstruct(nstruc,losse,losst,N,direction)
% Function to find the best structure for all model orders
% This function is created for the course Stochastic Adaptive Control, then
% modified for the master thesis

% Estimation data
disp('Estimations data')
beststre=[]; bestlosse=[];
np=sum(nstruc(:,1:3)')';
figure(40)
plot(np,losse,'+'); grid; title('Loss - Estimation set')
% ylim([1 1.5])
xlabel('Number of parameters');
% The next iteration is for finding the smallest loss structure
%   variation for all the given structure sizes
for i=min(np):max(np),
    idx=find(i==np);
    structures=nstruc(idx,:);
    losses=losse(idx);
    [loss,index]=min(losses);
    bestlosse=[bestlosse; loss];
    beststre=[beststre; structures(index,:)];
    %disp([loss structures(index,:)])
end
np=min(np):max(np); np=np(:);
hold on
plot(np,bestlosse,'-');
hold off

% ---- Test data -----
disp('Test data')
disp('    Loss               structure ')
bestlosst=[]; beststrt = [];
np=sum(nstruc(:,1:3)')';
figure(40)
plot(np,losst,'+'); grid; title('Loss - Test data');
% ylim([1 1.15])
xlabel('# of parms');
for i=min(np):max(np),
    idx=find(i==np);
    structures=nstruc(idx,:);
    losses=losst(idx);
    [loss,index]=min(losses);
    bestlosst=[bestlosst; loss];
    beststrt=[beststrt; structures(index,:)];
    %disp([loss structures(index,:)]);
end
np=min(np):max(np); np=np(:);
hold on
plot(np,bestlosst,'-');
hold off
loss=[bestlosse bestlosst np]; % and ns contain information to be used later

% Parameters
disp('      Diff                Loss               #parm')
disp([[NaN NaN; diff(loss(:,1:2))] loss])
h = figure(40)
hold on
plot(np(1:end-1),loss(1:end-1,1:end-1),'*-');
hold off
if direction == 2
    legend({'Estimation data, forward','Test data, forward',...
        'Estimation data, backward','Test data, backward'},'Location','east');
%     ylim([0.95 1.15]), grid;
    title('Loss function'); xlabel('Number of parameters');
    % saveas(h,'lossfunctions.png')
end


% -- Ftest --
n=length(loss);
ftable=zeros(n-1,n-1); 
lossA=loss(1:end-1,1); npA=loss(1:end-1,end);
lossB=loss(2:end,1);   npB=loss(2:end,end);
for i=1:n-1,
    for j=i:n-1,
        % to remember: N = length(Ze)
        dloss=(lossA(i)-lossB(j))/lossB(j);
        dpar=(N-npB(j))/(npB(j)-npA(i));
        ftest=dloss*dpar;
        if ftest>0,
            ftable(i,j)=round(fcdf(ftest,npB(j)-npA(i),N-npB(j))*100);
        end
    end;
end;
disp('Models'); disp(' ');
disp([np beststrt]);
disp(' ')
disp('F test');
disp(' ')
disp([ NaN np(2:end)'; np(1:end-1) ftable]);
disp(' ')

% -- Aic, Bic FPE --
np=loss(:,end); loss=loss(:,1);
ic=zeros(n,3);
for i=1:n,
    ic(i,1)=(1+2*np(i)/N)*loss(i);       % AIC
    ic(i,2)=(1+log(N)*np(i)/N)*loss(i);  % BIC
    ic(i,3)=(N+np(i))*loss(i)/(N-np(i)); % FPE
end
figurehandle = figure(30);
subplot(131)
% h5=figure(5);
grid;
if direction == 1
    hold on
    plot(np(1:end-1),ic(1:end-1,1),'o-');
    plot(np(1:end-1),ic(1:end-1,2),'o-');
    plot(np(1:end-1),ic(1:end-1,3),'o--');
    hold off
else
    hold on
    plot(np(1:end-1),ic(1:end-1,1),'*-');
    plot(np(1:end-1),ic(1:end-1,2),'*-');
    plot(np(1:end-1),ic(1:end-1,3),'*--');
    hold off
    title('AIC, BIC, FPE'); xlabel('Number of parameters');
    legend({'AIC forward','BIC forward','FPE forward',...
        'AIC backward','BIC backward','FPE backward'},'Location','east')
%     ylim([1 1.15])
    grid;
    % saveas(h5,'infcriteria.png')
end

% Display of information criteria
[mi,im]=min(ic);
disp(' '); disp('AIC BIC FPE'); disp(' ');
disp(np(im)')
AIC = np(im(1));
BIC = np(im(2));
FPE = np(im(3));




