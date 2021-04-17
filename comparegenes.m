mg = microgrid(2,3600);
N = 24;
uopt = 4;
chromosomes = 3;
genlength = 11;

[x0,cbuy,csale,cprod,Pload,Pres] = mg.representativeSample(N);
z = {x0,cbuy,csale,cprod,Pload,Pres};

controller1 = microgridMPCController(mg,N);
u1 = controller1.control(z);
utilde = u1{1}(uopt:4:end);
utilde1 = u1{1}(1:4:end);
utilde2 = u1{1}(2:4:end);
utilde3 = u1{1}(3:4:end);

constants = [mg.xlb,mg.xub,mg.ulb(1)',mg.ulb(3:4)',mg.uub(1)',mg.uub(3:4)',0.5,1,2,-1];
v(1,1) = x0;
for k = 1:N
    v(2,k) = cbuy(k);
    v(3,k) = csale(k);
    v(4,k) = cprod(k);
    v(5,k) = Pload(k);
    v(6,k) = Pres(k);
    if k > 1
        [M1,M2,M3] = mg.NstepDynamicsMatrices(k-1);
        xt = M1*[u1{1}(1:(k-1)*4);u1{2}(1:(k-1)*4);u1{3}(1:(k-1)*2)] + M2*x0 + M3;
        v(1,k) = xt;
        v(7,k) = cbuy(k-1);
        v(8,k) = csale(k-1);
        v(9,k) = cprod(k-1);
        v(10,k) = Pload(k-1);
        v(11,k) = Pres(k-1);
    end
end

for m = 0:chromosomes-1
    chrom = gen1((m*genlength)*3+1:3*(m+1)*genlength);
    expressions{m+1} = Tree.fromGen(char(chrom));
end
f = zeros(N,3);
for m = 1:3
    for k = 1:N
        f(k,m) = expressions{m}.evaluate(constants,v(:,k));
    end
end
uhat = f*(f\utilde);

for m = 0:chromosomes-1
    chrom = gen2((m*genlength)*3+1:3*(m+1)*genlength);
    expressions{m+1} = Tree.fromGen(char(chrom));
end
f2 = zeros(N,3);
for m = 1:3
    for k = 1:N
        f2(k,m) = expressions{m}.evaluate(constants,v(:,k));
    end
end
uhat2 = f2*(f2\utilde);

figure
plot(1:N,utilde,1:N,uhat,'r--',1:N,uhat2,'k--')
title('Input parametrization')
legend('optimal input','gen1','gen2')
ylabel('power exchanged(kW)')
xlabel('k')

annih1 = eye(N)-f*pinv(f);
se1 = utilde'*annih1*utilde

annih2 = eye(N)-f2*pinv(f2);
se2 = utilde'*annih2*utilde