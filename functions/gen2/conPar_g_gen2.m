function g = conPar_g_gen2(x,Prob)
global US_A

delta = Prob.user.delta;
cs = Prob.user.cs;
N = Prob.user.N;

zdM = [1 0 0 0;0 0 0 1];
zuM = [1 0 0 0;0 0 0 1];

dU = zeros(N*4,9);

%du1 = zeros(1,9);
du2 = zeros(1,9);
du3 = zeros(1,9);
du4 = zeros(1,9);


for i = 1:N
    du4(1:3) = [US_A.f1(i),US_A.f2(i),US_A.f3(i)];
    
    du2(4:6) = [US_A.f4(i),US_A.f5(i),US_A.f6(i)];

    du3(7:9) = [US_A.f4(i),US_A.f5(i),US_A.f6(i)];

    du1 = du4+du2+du3;

    dU((i-1)*4+1:4*i,:) = [du1;du2;du3;du4]; 
end

dZ = zeros(9,1);
for i = 1:9
    Z = (kron(eye(N),zuM)*dU(:,i)).*(kron(eye(N),zdM)*delta);
    dZ(i) = cs(end-2*N+1:end)'*Z;
end

g = (cs(1:4*N)'*dU)'+dZ;
end

