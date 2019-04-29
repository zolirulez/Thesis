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
    eval([constantfields{i} '=constants.' constantfields{i} ]);
end

if ~all(u==0)
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
end
if ~all(x==0)
    % State values
    DVAValue = x(1);
    p1Value = x(2); % p1Value = initx(2);
    h1Value = x(3);
    d1Value = initx(4);
    TA2Value = x(5); % TA2Value = initx(5);
    Dm21Value = x(6);
    p2Value = x(7);
    h2Value = x(8);
    d2Value = initx(9);
    TA1Value = x(10);
    BPValue = x(11); % BPValue = initx(11);
    DmVValue = x(12);
    pRValue = x(13);
    hRValue = x(14);
    dRValue = initx(15);
    DmGValue = x(16);
    delta_h2Value = initx(17);
    DmQValue = x(18);
else
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

% Substitutions
ABCD = eval(['ABCD(' value ')']);
A = ABCD(1:nx,1:nx);
B = ABCD(1:nx,nx+1:nx+nu);
C = ABCD(nx+1:nx+ny,1:nx);
D = ABCD(nx+1:nx+ny,nx+1:nx+nu);
ABCD = c2d(ss(A,B,C,D),1);

ABCD = [reshape(ABCD.A,nx*nx,1); ...
    reshape(ABCD.B,nx*nu,1); ...
    reshape(ABCD.C,ny*nx,1); ...
    reshape(ABCD.D,ny*nu,1)];