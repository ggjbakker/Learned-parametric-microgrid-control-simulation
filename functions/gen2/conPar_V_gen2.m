function V = conPar_V_gen2(x,Prob)
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

zdM = [1 0 0 0;0 0 0 1];
zuM = [1 0 0 0;0 0 0 1];

f1 = zeros(N,1);
f2 = zeros(N,1);
f3 = zeros(N,1);
f4 = zeros(N,1);
f5 = zeros(N,1);
f6 = zeros(N,1);
f7 = zeros(N,1);
f8 = zeros(N,1);
f9 = zeros(N,1);

f1(1) = ulb(3)-csale(1)+2*x0;
f2(1) = csale0*x0*cbuy0^2*(cbuy(1)-1);
f3(1) = csale(1)^2*x0*(1+uub(4));

f4(1) = csale(1);
f5(1) = csale(1)*cprod(1)-2*cbuy(1);
f6(1) = cprod(1);

% f7(1) = cprod(1);
% f8(1) = csale(1);
% f9(1) = (Pload0+Pres0+u0(2)+u0(3))/1500;

xs = x0;
u4 = x(1)*f1(1) + x(2)*f2(1) + x(3) * f3(1);
u2 = f4(1)*x(4) + x(5)*f5(1) + x(6)*f6(1);
u3 = f4(1)*x(7) + x(8)*f5(1) + x(9)*f6(1);
%u4 = f7(1)*x(7) + f8(1)*x(8) + x(9)*f9(1)+x(10);

u1 = u4+u2+u3+Pres(1)-Pload(1);

U = [u1;u2;u2;u4;zeros((N-1)*4,1)];

for i = 1:N-1
    ru = 4*(i-1)+1:4*(i-1)+4;
    rd = 4*(i-1)+1:4*(i-1)+4;
    xs = mA*xs + mB1 * U(ru) + mB2*delta(rd) + mB3*((zdM*delta(rd)).*(zuM*U(ru)));
    
    f1(i+1) = ulb(3)-csale(i+1)+2*xs;
    f2(i+1) = csale(i)*xs*cbuy(i)^2*(cbuy(i+1)-1);
    f3(i+1) = csale(i+1)^2*xs*(1+uub(4));

    f4(i+1) =  csale(i+1);
    f5(i+1) = csale(i+1)*cprod(i+1)-2*cbuy(i+1);
    f6(i+1) = cprod(i+1);

%     f7(i+1) = cprod(i+1);
%     f8(i+1) = csale(i+1);
%     f9(i+1) = (Pload(i)+Pres(i)+u2+u3)/1500;
    
    u4 = x(1)*f1(i+1) + x(2)*f2(i+1) + x(3) * f3(i+1);
    u2 = f4(i+1)*x(4) + x(5)*f5(i+1) + x(6)*f6(i+1);
    u3 = f4(i+1)*x(7) + x(8)*f5(i+1) + x(9)*f6(i+1);
%     u4 = f7(i+1)*x(7) + f8(i+1)*x(8) + x(9)*f9(i+1)+x(10);
    
    u1 = u4+u2+u3+Pres(i+1)-Pload(i+1);

    U(ru+4) = [u1;u2;u3;u4];
end

Z = (kron(eye(N),zuM)*U).*(kron(eye(N),zdM)*delta);

V = {};
V{1} = U;
V{2} = delta;
V{3} = Z;
end

