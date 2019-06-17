TA0 = U(1:7000,12);
TBP = Y(1:7000,2);
[cc, lags] = xcov(TA0,TBP);
figure, plot(lags,cc)