function ABCDQ = LTVsystemDescription(u,x,y,w)

% Use x or y (x means feedback, y means rough assumptions)
% feedback = 0;

constants = load('constants.mat');
constantfields = fieldnames(constants);
for i=1:length(constantfields)
    eval([constantfields{i} '=constants.' constantfields{i} ';']);
end

% Input values
CRAValue = u(1);
BPValue = u(2);
CRVValue = u(3);
CRITValue = u(4);
delta_hHRValue = u(5);
dBPValue = u(6);
dGValue = u(7);
hGValue = u(8);
hLValue = u(9);
TA0Value = u(10);
hHRValue = u(11);
DmQValue = u(12);

% States MODIFIED (TODO)
p1Value = y(1);
h1Value = x(2);
d1Value = x(3);
TA1Value = x(4);
pRValue = y(3);
hRValue = y(4);
dRValue = x(7);
delta_hValue = x(8);
% Tables
p1idx = max(1,[find(p1Value < Pbig)-1 51]);
pRidx = max(1,[find(pRValue < Pbig)-1 51]);
h1idx = max(1,[find(h1Value < Hbig)-1 51]);
hRidx = max(1,[find(hRValue < Hbig)-1 51]);
delta_ph1Value = paramvectorP(1,p1idx(1),h1idx(1));
delta_phRValue = paramvectorP(1,pRidx(1),hRidx(1));
delta_pd1Value = paramvectorP(2,p1idx(1),h1idx(1));
delta_pdRValue = paramvectorP(2,pRidx(1),hRidx(1));
delta_Th1Value = paramvectorT(1,p1idx(1),h1idx(1));
delta_Td1Value = paramvectorT(2,p1idx(1),h1idx(1));

% save('deltaValues','delta_Th1Value','delta_Td1Value');

% Parameter estimation results
s0Value = w(1);
kValue = w(2);

% Substitutions
value = [strrep('BP,CRA,CRIT,CRV,DmQ,TA0,TA1,d1,dBP,dG,dR,delta_h,delta_Td1,delta_Th1,delta_pd1,delta_ph1,delta_pdR,delta_phR,h1,hG,hHR,hL,k,p1,pR,s0',',','Value,') 'Value'];
ABCD = eval(['ABCD(' value ');']);
A = ABCD(1:nx,1:nx);
B = ABCD(1:nx,nx+1:nx+nu);
C = ABCD(nx+1:nx+ny,1:nx);
D = ABCD(nx+1:nx+ny,nx+1:nx+nu);
Ts = 1;
[A,B,Q] = c2dn(A,B,Lq,Ts);

ABCDQ = [reshape(A,nx*nx,1); ...
    reshape(B,nx*nu,1); ...
    reshape(C,ny*nx,1); ...
    reshape(D,ny*nu,1); ...
    reshape(Q,nx*nx,1)];