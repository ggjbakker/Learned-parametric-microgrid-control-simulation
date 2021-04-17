
dt = 120;
dtC = 1800;
N = 48;
x0 = 110;
mg = microgrid(2,dtC,0.9,0.9,500,5,-150,150,250,20,-1000,1000,x0);
uopt = 4;
chromosomes = 3;
genlength = length(gen)/3/chromosomes;

[x0,cbuy,csale,cprod,Pload,Pres] = mg.representativeSample(N);
z = {x0,cbuy,csale,cprod,Pload,Pres};

controller1 = microgridMPCController(mg,N);
[~, result] = controller1.control(z);
u1 = result.V;
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
    chrom = gen((m*genlength)*3+1:3*(m+1)*genlength);
    expressions{m+1} = Tree.fromGen(char(chrom));
end
f = zeros(N,3);
for m = 1:3
    for k = 1:N
        f(k,m) = expressions{m}.evaluate(constants,v(:,k));
    end
end

uhat = f*(f\utilde);
f2 = zeros(N,3);
for k = 1:N
    f2(k,1) = -cprod(k);
    f2(k,2) = -csale(k);
    if k > 1
        f2(k,3) = -Pload(k-1) + Pres(k-1)+utilde2(k-1)+utilde3(k-1);
    end
    
end
uhat2 = f2*(f2\utilde);
figgi = figure;
figgi.Position =  [100 100 800 400];
plot(1:N-1,utilde(2:N),1:N-1,uhat(2:N),'r-.',1:N-1,uhat2(2:N),'k--')
title({'Example parametrization $P^*_\mathrm{grid}$'},'Interpreter','latex','FontSize',16)
legend({'$P^*_\mathrm{grid}$','Expression tree','Hand crafted'},'Interpreter','latex','FontSize',10)
ylabel({'$P_\mathrm{grid}(\mathrm{kW})$'},'Interpreter','latex','FontSize',14)
xlabel('$h$','Interpreter','latex','FontSize',14)

annih1 = eye(N)-f*pinv(f);
fitnesslearned = utilde'*annih1*utilde/(N-3)

annih2 = eye(N)-f2*pinv(f2);
fitnesshandcrafted = utilde'*annih2*utilde/(N-3)