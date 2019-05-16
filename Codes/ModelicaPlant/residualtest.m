load longsimulation_faulty
delta_h = record(8,:);
% Note: dont start it before 350, because of the delay operations
start = 2001;
resrecord = NaN(2,finish-start+1);
W = [5000; 6000; 1e3; 1e4; 1e3]*[1 10]; 
P = diag(max(W')')/10; 
for it = start:finish
    TA0 = U(it,10);
    THR = CoolProp.PropsSI('T','H',U(it,11),'P',Y(it,1),'CO2');
    CRIT = U(it-1,4)*VValues(3)*U(it-1,7);
    DV = U(it-1,1)*DVValues(2);
    DmV = U(it,3)*KvValues*sqrt(U(it,6)*(Y(it,1) - Y(it,3)));
    DQ = DmV*(U(it,11) - hBP);
    DT = TBP - TA0;
    sigma = DQ/DT;
    phi = [1; DV; TA0-273.15; CRIT; THR-273.15];
    res = [sigma delta_h(it-start+1)] - phi'*W;
    schur = 1 + phi'*P*phi;
    K = P*phi/schur ;
    W = W + K*res; 
    P = P - K*schur*K';
    resrecord(:,it-start+1) = res;
end
figure(5)
plot(start:finish,resrecord')