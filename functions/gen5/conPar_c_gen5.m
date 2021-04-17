function c = conPar_c_gen5(x,Prob)
global US_B
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
csale0 = Prob.user.csale0;
cbuy0 = Prob.user.cbuy0;

mA = Prob.user.mA;
mB1 = Prob.user.mB1;
mB2 = Prob.user.mB2;
mB3 = Prob.user.mB3;

cs = Prob.user.cs;

xub = Prob.user.xub;
xlb = Prob.user.xlb;

uub = Prob.user.uub;
ulb = Prob.user.ulb;

N = Prob.user.N;

A1 = Prob.user.A1;
A2 = Prob.user.A2;
A3 = Prob.user.A3;
A4 = Prob.user.A4;

zdM = [1 0 0 0;0 0 0 1];
zuM = [1 0 0 0;0 0 0 1];

US_B = struct;

US_B.f1 = zeros(N,1);
US_B.f2 = zeros(N,1);
US_B.f3 = zeros(N,1);
US_B.f4 = zeros(N,1);
US_B.f5 = zeros(N,1);
US_B.f6 = zeros(N,1);

US_B.f1(1) = csale(1)-cprod(1);
US_B.f2(1) = uub(3)*(csale(1)*cbuy0-0.5)*(cbuy(1)-1);
US_B.f3(1) = (2*cbuy0+Pres(1))*(cbuy0+x0+Pres(1)+uub(4));

US_B.f4(1) = csale(1);
US_B.f5(1) = cprod(1)^2-cprod(1)*ulb(1);
US_B.f6(1) = cprod(1);

xs = x0;
u4 = x(1)*US_B.f1(1) + x(2)*US_B.f2(1) + x(3) * US_B.f3(1);
u2 = US_B.f4(1)*x(4) + x(5)*US_B.f5(1) + x(6)*US_B.f6(1);
u3 = US_B.f4(1)*x(7) + x(8)*US_B.f5(1) + x(9)*US_B.f6(1);

u1 = u4+u2+u3+Pres(1)-Pload(1);

U = [u1;u2;u2;u4;zeros((N-1)*4,1)];

for i = 1:N-1
    ru = 4*(i-1)+1:4*(i-1)+4;
    rd = 4*(i-1)+1:4*(i-1)+4;
    xs = mA*xs + mB1 * U(ru) + mB2*delta(rd) + mB3*((zdM*delta(rd)).*(zuM*U(ru)));
    
    US_B.f1(i+1) = csale(i+1)-cprod(i+1);
    US_B.f2(i+1) = uub(3)*(csale(i+1)*cbuy(i)-0.5)*(cbuy(i+1)-1);
    US_B.f3(i+1) = (2*cbuy(i)+Pres(i+1))*(cbuy(i)+xs+Pres(i+1)+uub(4));

    US_B.f4(i+1) = csale(i+1);
    US_B.f5(i+1) = cprod(i+1)^2-cprod(i+1)*ulb(1);
    US_B.f6(i+1) = cprod(i+1);

    
    u4 = x(1)*US_B.f1(i+1) + x(2)*US_B.f2(i+1) + x(3) * US_B.f3(i+1);
    u2 = US_B.f4(i+1)*x(4) + x(5)*US_B.f5(i+1) + x(6)*US_B.f6(i+1);
    u3 = US_B.f4(i+1)*x(7) + x(8)*US_B.f5(i+1) + x(9)*US_B.f6(i+1);
    
    u1 = u4+u2+u3+Pres(i+1)-Pload(i+1);

    U(ru+4) = [u1;u2;u3;u4];
end

Z = (kron(eye(N),zuM)*U).*(kron(eye(N),zdM)*delta);
c = [A1 * [U;delta;Z];A2 * [U;delta;Z];A3 * [U;delta;Z];A4 * [U;delta;Z]];
end

