function dx = dynamicsf(t,x,n,Design_Parametersn,B,C,D,G,J,fx,fy)

q = x(1:n);
dq = x(n+1:end);
%#codegen
tx = transpose(x);
M1 = D(Design_Parametersn,tx);
C1 = C(Design_Parametersn,tx);
G1 = G(Design_Parametersn,tx);
B1 = B(Design_Parametersn,tx);
%u = zeros(n,1);
u = transpose(J(Design_Parametersn,tx))*[fx ;fy];
%u(3) = 0.1;
B0 = 100*eye(n,n);%friction
ddq = M1\(-C1*dq-B0*dq-G1+B1*u);%
dx = [dq;ddq];
end