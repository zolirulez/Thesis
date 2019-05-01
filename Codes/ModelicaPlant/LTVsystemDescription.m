function ABCD = LTVsystemDescription(u,x,y,feedback)

% Use x or y (x means feedback, y means rough assumptions)
% feedback = 0;

constants = load('constants.mat');
constantfields = fieldnames(constants);
for i=1:length(constantfields)
    eval([constantfields{i} '=constants.' constantfields{i} ';']);
end

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


if ~feedback
    p2Value = y(1);
    hBPValue = y(2);
    pRValue = y(3);
    hRValue = y(4);
    hHRValue = y(5);
    TA1Value = y(6);
    DmQValue = y(7);
    
    % Used state values
    DVAValue = c.MxDVAValue*CRAValue;
    p1Value = p2Value + 1.5e5;
    h1Value = (hBPValue + initx(17) + hHRValue)/2;
    d1Value = initx(4);
    sValue = c.s0Value + c.kValue*DVAValue;
    wValue = c.dAValue*DVAValue*c.cpValue/sValue;
    h2Value = hBPValue + initx(17);
    
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
    % State values
    TA2Value = 1/(wValue+1)*(delta_Th1Value*h1Value + delta_Td1Value*d1Value) + 1/(1/wValue+1)*TA1Value;
    Dm21Value = 1/c.RValue*sqrt(d1Value*(p1Value-p2Value));
    %     p2Value = initx(7);
    d2Value = initx(9);
    %     TA1Value = initx(10);
    BPValue = BPRValue;
    DmVValue = CRVValue*c.KvVValue*sqrt(dBPValue*(p2Value - pRValue));
    %     pRValue = initx(13);
    %     hRValue = initx(14);
    dRValue = initx(15);
    DmGValue = dGValue*c.VGValue*CRITValue*c.MxfITValue;
    delta_h2Value = initx(17);
    %     DmQValue = initx(18);
else
    % States
    DVAValue = x(1);
    p1Value = x(2);
    h1Value = x(3);
    d1Value = x(4);
    TA2Value = x(5);
    Dm21Value = x(6);
    p2Value = x(7);
    h2Value = x(8);
    d2Value = x(9);
    TA1Value = x(10);
    BPValue = x(11);
    DmVValue = x(12);
    pRValue = x(13);
    hRValue = x(14);
    dRValue = x(15);
    DmGValue = x(16);
    delta_h2Value = x(17);
    DmQValue = x(18);
    % Tables
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