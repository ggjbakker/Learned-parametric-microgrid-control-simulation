function f = conPar_f(x,Prob)
global US_A
u0 = Prob.user.u0;
x0 = Prob.user.x0;
delta = Prob.user.delta;

cprod = Prob.user.cprod;
csale = Prob.user.csale;
Pload = Prob.user.Pload;
Pres = Prob.user.Pres;
cbuy = Prob.user.cbuy;
Pload0 = Prob.user.Pload0;
Pres0 = Prob.user.Pres0;

mA = Prob.user.mA;
mB1 = Prob.user.mB1;
mB2 = Prob.user.mB2;
mB3 = Prob.user.mB3;

cs = Prob.user.cs;

xub = Prob.user.xub;
xlb = Prob.user.xlb;

N = Prob.user.N;

zdM = [1 0 0 0;0 0 0 1];
zuM = [1 0 0 0;0 0 0 1];

US_A = struct;

US_A.f1 = zeros(N,1);
US_A.f2 = zeros(N,1);
US_A.f3 = zeros(N,1);
US_A.f4 = zeros(N,1);
US_A.f5 = zeros(N,1);
US_A.f6 = zeros(N,1);
US_A.f7 = zeros(N,1);
US_A.f8 = zeros(N,1);
US_A.f9 = zeros(N,1);

US_A.f1(1) = cbuy(1);
US_A.f2(1) = (Pload0+Pres0+u0(2)+u0(3))/1500;
US_A.f3(1) = (0.5*(xub-xlb) - x0)/125;

US_A.f4(1) = Pload(1)/750;
US_A.f5(1) = cbuy(1);
US_A.f6(1) = (xub - x0)/250;

xs = x0;

u2 = US_A.f4(1)*x(1) + x(2)*US_A.f5(1) + x(3)*US_A.f6(1);
u3 = US_A.f4(1)*x(4) + x(5)*US_A.f5(1) + x(6)*US_A.f6(1);
u4 = US_A.f1(1)*x(7) + US_A.f2(1)*x(8) + x(9)*US_A.f3(1)+x(10);

u1 = u4+u2+u3+Pres(1)-Pload(1);

U = [u1;u2;u2;u4;zeros((N-1)*4,1)];

for i = 1:N-1
    ru = 4*(i-1)+1:4*(i-1)+4;
    
    US_A.f1(i+1) = cbuy(i+1);
    US_A.f2(i+1) = (Pload(i)+Pres(i)+u2+u3)/1500;
    US_A.f3(i+1) = (0.5*(xub-xlb) - xs)/125;

    US_A.f4(i+1) = Pload(i+1)/750;
    US_A.f5(i+1) = cbuy(i+1);
    US_A.f6(i+1) = (xub - xs)/250;

    u2 = US_A.f4(i+1)*x(1) + x(2)*US_A.f5(i+1) + x(3)*US_A.f6(i+1);
    u3 = US_A.f4(i+1)*x(4) + x(5)*US_A.f5(i+1) + x(6)*US_A.f6(i+1);
    u4 = US_A.f1(i+1)*x(7) + US_A.f2(i+1)*x(8) + x(9)*US_A.f3(i+1)+x(10);
    
    u1 = u4+u2+u3+Pres(i+1)-Pload(i+1);

    U(ru+4) = [u1;u2;u3;u4];
end

Z = (kron(eye(N),zuM)*U).*(kron(eye(N),zdM)*delta);

f = cs'*[U;delta;Z];
end

