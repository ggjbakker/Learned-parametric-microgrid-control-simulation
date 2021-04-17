addpath('.\classes')
addpath('.\functions')

dtC = 1800;
x0 = 110;
mgrid = microgrid(2,dtC,0.9,0.9,495,5,-145,145,245,21,-1000,1000,x0);
N = 12;
Nz = 2000;
Nztest = 20000;
Noptsteps = 50;
nv = 11;
uopt = 4;
%%S
controller = microgridMPCController(mgrid,N);
xs = mgrid.representativeSet(N,Nztest);
dataset = cell(N,10);

f = waitbar(0,'creating test dataset');
vtest = zeros(nv,N,Nztest);
v = zeros(nv,N,Nz);
u1 = zeros(N,Nz);
u2 = zeros(N,Nz);
u4 = zeros(N,Nz);
u1test = zeros(N,Nz);
u2test = zeros(N,Nz);
u4test = zeros(N,Nz);
data1 = {};
data2 = {};
data4 = {};
testdata1 = {};
testdata2 = {};
testdata4 = {};
for i = 1:size(xs,1)
    [~,result] = controller.control(xs(i,:));
    solution = result.V;
    dataset(i,1:6) = xs(i,:);
    dataset(i,7:9) = solution;
    vtest(1,1,i) = cell2mat(xs(i,1));
    u1test(:,i) = solution{1}(1:4:end);
    u2test(:,i) = solution{1}(2:4:end);
    u4test(:,i) = solution{1}(4:4:end);
    for k = 2:N
        [M1,M2,M3] = mgrid.NstepDynamicsMatrices(k-1);
        x = M1*[solution{1}(1:(k-1)*4);solution{2}(1:(k-1)*4);solution{3}(1:(k-1)*2)] + M2*cell2mat(xs(i,1)) + M3;
        vtest(1,k,i) = x;
    end
    for k = 2:6
        vtest(k,:,i) = cell2mat(xs(i,k))';
    end
    for k = 7:11
        pv = cell2mat(xs(i,k-5))';
        vtest(k,:,i) = [0 pv(1:end-1)];
    end
    waitbar(i/size(xs,1),f,'creating test dataset')
end
close(f)
dataset = cellfun(@transpose,dataset,'UniformOutput',false);

c = [mgrid.xlb,mgrid.xub,mgrid.ulb(1)',mgrid.ulb(3:4)',mgrid.uub(1)',mgrid.uub(3:4)',0.5,1,2,-1];

testdata1.v = vtest;
testdata1.u = u1test;
pool1 = Pool(20,7,3,testdata1,c);
pool1.populate();

testdata2.v = vtest;
testdata2.u = u2test;
pool2 = Pool(20,7,3,testdata2,c);
pool2.populate();

testdata4.v = vtest;
testdata4.u = u4test;
pool4 = Pool(20,7,3,testdata4,c);
pool4.populate();

j = 1;
k = 0;
%%

while j < 100000
    if mod(j+2,3) == 0
        k = k+1;
        pool1.nsamples = size(testdata1.v,3);
        pool1.data = testdata1;
        fitness = pool1.evaluatePopulation();
        meanfitness1(k) = mean(fitness(:,2));
        best1(k) = min(fitness(:,2));
        fprintf('\nTest result u1 for iteration %i \n',j);
        fprintf('Mean u1: %d \n',meanfitness1(k));
        fprintf('Best u1: %d \n \n',best1(k));

        pool2.nsamples = size(testdata2.v,3);
        pool2.data = testdata2;
        fitness = pool2.evaluatePopulation();
        meanfitness2(k) = mean(fitness(:,2));
        best2(k) = min(fitness(:,2));
        fprintf('Mean u2: %d \n',meanfitness2(k));
        fprintf('Best u2: %d \n \n',best2(k));

        pool4.nsamples = size(testdata4.v,3);
        pool4.data = testdata4;
        fitness = pool4.evaluatePopulation();
        meanfitness4(k) = mean(fitness(:,2));
        best4(k) = min(fitness(:,2));
        fprintf('Mean u4: %d \n',meanfitness4(k));
        fprintf('Best u4: %d \n \n',best4(k));

        if best1(k) == min(best1)
            save('.\datasets\Genpar5r2\bestpool1.mat','pool1')
        end
        if best2(k) == min(best2)
            save('.\datasets\Genpar5r2\bestpool2.mat','pool2')
        end
        if best4(k) == min(best4)
            save('.\datasets\Genpar5r2\bestpool4.mat','pool4')
        end

        save('.\datasets\Genpar5r2\checkpoint.mat','meanfitness1','best1','pool1',...
            'meanfitness2','best2','pool2', ...
            'meanfitness4','best4','pool4', ...
            'testdata1','testdata2','testdata4','j','k')
    end
    
    xs = mgrid.representativeSet(N,Nz);
    for i = 1:size(xs,1)
        [~,result] = controller.control(xs(i,:));
        solution = result.V;
        dataset(i,1:6) = xs(i,:);
        dataset(i,7:9) = solution;
        v(1,1,i) = cell2mat(xs(i,1));
        u1(:,i) = solution{1}(1:4:end);
        u2(:,i) = solution{1}(2:4:end);
        u4(:,i) = solution{1}(4:4:end);
        for m = 2:N
            [M1,M2,M3] = mgrid.NstepDynamicsMatrices(m-1);
            x = M1*[solution{1}(1:(m-1)*4);solution{2}(1:(m-1)*4);solution{3}(1:(m-1)*2)] + M2*cell2mat(xs(i,1)) + M3;
            v(1,m,i) = x;
        end
        for m = 2:6
            v(m,:,i) = cell2mat(xs(i,m))';
        end
        for m = 7:11
            pv = cell2mat(xs(i,m-5))';
            v(m,:,i) = [0 pv(1:end-1)];
        end
    end
    data1.v = v;
    data1.u = u1;
    pool1.data = data1;
    pool1.nsamples = size(data1.v,3);
    
    data2.v = v;
    data2.u = u2;
    pool2.data = data2;
    pool2.nsamples = size(data2.v,3);

    data4.v = v;
    data4.u = u4;
    pool4.data = data4;
    pool4.nsamples = size(data4.v,3);
    
    fprintf('Evolving genes step: %4d\n',1);
    for i = 1:Noptsteps
        fprintf(1,'\b\b\b\b\b%4d\n',i)
        [m1,b1] = pool1.step();
        [m2,b2] = pool2.step();
        [m4,b4] = pool4.step();
    end
    j = j+1;
end
%%
xs = mgrid.representativeSet(N,1000);
f = waitbar(0,'creating dataset');
for i = 1:size(xs,1)
    [~,result] = controller.control(xs(i,:));
    solution = result.V;
    dataset(i,1:6) = xs(i,:);
    dataset(i,7:9) = solution;
    v(1,1,i) = cell2mat(xs(i,1));
    u1(:,i) = solution{1}(1:4:end);
    u2(:,i) = solution{1}(2:4:end);
    u4(:,i) = solution{1}(4:4:end);
    for k = 2:N
        [M1,M2,M3] = mgrid.NstepDynamicsMatrices(k-1);
        x = M1*[solution{1}(1:(k-1)*4);solution{2}(1:(k-1)*4);solution{3}(1:(k-1)*2)] + M2*cell2mat(xs(i,1)) + M3;
        v(1,k,i) = x;
    end
    for k = 2:6
        v(k,:,i) = cell2mat(xs(i,k))';
    end
    for k = 7:11
        pv = cell2mat(xs(i,k-5))';
        v(k,:,i) = [0 pv(1:end-1)];
    end
    waitbar(i/size(xs,1),f,'creating dataset')
end
close(f)
data.v = v;
data.u = u4;
pool4.nsamples = size(data.v,3);
pool4.data = data;
fitness = pool4.evaluatePopulation();
meanfitness = mean(fitness(:,2))
best = min(fitness(:,2))