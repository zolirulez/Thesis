function ABCD = LTVsystemDescription(u,x)

% p2Value = y(1);
% hBPValue = y(2);
% pRValue = y(3);
% hRValue = y(4);
% hHRValue = y(5);
% % TA1Value = y(6);
% DmQValue = y(7);

constants = load('constants.mat')
constantfields = fieldnames(constants);
for i=1:length(constantfields)
    eval([constantfields{i} '=constants.' constantfields{i} ';']);
end

if 1%~all(u==0)
    % Input values
    CRAValue = u(1);
    BPRValue = u(2);
    CRVValue = u(3);
    CRITValue = u(4);
    delta_hHRValue = u(5);
    dBPValue = u(6);
    dGValue = u(7);
    hGValue = u(8);
    hLValue = u(9);
    TA0Value = u(10);
    hMTValue = u(11);
    % State values
    DVAValue = initx(1);
    p1Value = initx(2); % p1Value = initx(2);
    h1Value = initx(3);
    d1Value = initx(4);
    TA2Value = initx(5); % TA2Value = initx(5);
    Dm21Value = initx(6);
    p2Value = initx(7);
    h2Value = initx(8);
    d2Value = initx(9);
    TA1Value = initx(10);
    BPValue = initx(11); % BPValue = initx(11);
    DmVValue = initx(12);
    pRValue = initx(13);
    hRValue = initx(14);
    dRValue = initx(15);
    DmGValue = initx(16);
    delta_h2Value = initx(17);
    DmQValue = initx(18);    
else
    % Input values
    CRAValue = initu(1);
    BPRValue = initu(2);
    CRVValue = initu(3);
    CRITValue = initu(4);
    delta_hHRValue = initu(5);
    dBPValue = initu(6);
    dGValue = initu(7);
    hGValue = initu(8);
    hLValue = initu(9);
    TA0Value = initu(10);
    hMTValue = initu(11);
    % State values
    DVAValue = initx(1);
    p1Value = initx(2); % p1Value = initx(2);
    h1Value = initx(3);
    d1Value = initx(4);
    TA2Value = initx(5); % TA2Value = initx(5);
    Dm21Value = initx(6);
    p2Value = initx(7);
    h2Value = initx(8);
    d2Value = initx(9);
    TA1Value = initx(10);
    BPValue = initx(11); % BPValue = initx(11);
    DmVValue = initx(12);
    pRValue = initx(13);
    hRValue = initx(14);
    dRValue = initx(15);
    DmGValue = initx(16);
    delta_h2Value = initx(17);
    DmQValue = initx(18);
end

% Table values
try
    p1idx = [find(p1Value < Pbig)-1 50];
    p2idx = [find(p2Value < Pbig)-1 50];
    pRidx = [find(pRValue < Pbig)-1 50];
    h1idx = [find(h1Value < Hbig)-1 50];
    h2idx = [find(h2Value < Hbig)-1 50];
    hRidx = [find(hRValue < Hbig)-1 50];
    delta_ph1Value = paramvectorP(1,p1idx(1),h1idx(1));
    delta_ph2Value = paramvectorP(1,p2idx(1),h2idx(1));
    delta_phRValue = paramvectorP(1,pRidx(1),hRidx(1));
    delta_pd1Value = paramvectorP(2,p1idx(1),h1idx(1));
    delta_pd2Value = paramvectorP(2,p2idx(1),h2idx(1));
    delta_pdRValue = paramvectorP(2,pRidx(1),hRidx(1));
    delta_Th1Value = paramvectorT(1,p1idx(1),h1idx(1));
    delta_Th2Value = paramvectorT(1,p2idx(1),h2idx(1));
    delta_Td1Value = paramvectorT(2,p1idx(1),h1idx(1));
    delta_Td2Value = paramvectorT(2,p2idx(1),h2idx(1));
catch
    disp('stop here')
end

% Substitutions
value = [strrep('BP,CRIT,CRV,DVA,Dm21,DmG,DmQ,DmV,TA0,TA1,TA2,d1,d2,dBP,dG,dR,delta_h2,delta_Td1,delta_Td2,delta_Th1,delta_Th2,delta_pd1,delta_pd2,delta_ph1,delta_ph2,delta_hHR,delta_pdR,delta_phR,h1,h2,hG,hL,hMT,p1,p2,pR',',','Value,') 'Value'];
ABCD = eval(['ABCD(' value ');']);
A = ABCD(1:nx,1:nx);
B = ABCD(1:nx,nx+1:nx+nu);
C = ABCD(nx+1:nx+ny,1:nx);
D = ABCD(nx+1:nx+ny,nx+1:nx+nu);
ABCD = c2d(ss(A,B,C,D),1);

ABCD = [reshape(ABCD.A,nx*nx,1); ...
    reshape(ABCD.B,nx*nu,1); ...
    reshape(ABCD.C,ny*nx,1); ...
    reshape(ABCD.D,ny*nu,1)];