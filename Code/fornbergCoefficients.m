ax=0;
bx=10;
m=999;
gridchoice='uniform';
x = xgrid(ax,bx,m,gridchoice);
#A1p = fdcoeffF(1, x(1), x(1:2));
#Aip =fdcoeffF(1, x(i), x((i-1):(i+1)));
#Aipp =fdcoeffF(2, x(i), x((i-1):(i+1)));

# 4th order accurate
A1p4 = fdcoeffF(1, x(1), x(1:5));
A2p4 = fdcoeffF(1, x(2), x(1:5));
A3p4 = fdcoeffF(1, x(3), x(1:5));
Aip4 = fdcoeffF(1, x(i), x((i-3):(i+1)));
Amp4 = fdcoeffF(1, x(m), x((m-5):m));

A1pp4 = fdcoeffF(2, x(1), x(1:5));
A2pp4 = fdcoeffF(2, x(2), x(1:5));
A3pp4 = fdcoeffF(2, x(3), x(1:5));
Aipp4 = fdcoeffF(2, x(i), x((i-2):(i+2)));
Am1pp4 = fdcoeffF(2, x(m-1), x((m-6):(m-1)));
Ampp4 = fdcoeffF(2, x(m), x((m-5):m));
