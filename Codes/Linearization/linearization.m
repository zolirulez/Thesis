% ------------------------ LINEARIZATION ----------------------------------                              
% Jacobians
A = jacobian(Dx,x);
B = jacobian(Dx,u);
GBC = jacobian(Dx,BC);
GdDer = jacobian(Dx,dDer);
GdDis = jacobian(Dx,dDis);
Gd = jacobian(Dx,d);
C = jacobian(y,x);
D = jacobian(y,u);
EBC = jacobian(y,BC);
EdDer = jacobian(y,dDer);
EdDis = jacobian(y,dDis);
Ed = jacobian(y,d);
mx4sub = struct('A',A,'B',B,'C',C,'Gd',Gd);
fields = fieldnames(mx4sub);

