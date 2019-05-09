record = [];
for it = 500:10000
    DmQ = U(it,3)*KvValues(1)*sqrt(U(it,6)*(Y(it,1) - Y(it,3))) -...
        U(it,7)*VValues(3)*U(it,4)*DVValues(2);
    record = [record DmQ];
end
figure(1)
plot(record(1,:))