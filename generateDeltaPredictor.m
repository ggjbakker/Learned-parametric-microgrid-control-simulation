addpath('.\classes')
addpath('.\functions')

dtC = 1800;
N = 48;
x0 = 110;
mgrid = microgrid(2,dtC,0.9,0.9,495,5,-145,145,245,21,-1000,1000,x0);
Nz = 20000;
numTrees = 10;
%%S
c = microgridMPCController(mgrid,N);
xs = mgrid.representativeSet(N,Nz);
dataset = cell(N,10);

f = waitbar(0,'creating dataset');
for i = 1:size(xs,1)
    [~,result] = c.control(xs(i,:));
    solution = result.V;
    dataset(i,1:6) = xs(i,:);
    dataset(i,7:9) = solution;
    waitbar(i/size(xs,1),f,'creating dataset')
end
dataset = cellfun(@transpose,dataset,'UniformOutput',false);
%%
binvars = mgrid.ndelta;
mdlsTree = {};

close(f)
f2 = waitbar(0,'Training tree');
for i = 1:N-1 
    
    response = cell(Nz,1);
    for k = 1:Nz
        response(k) = {num2str(logical(dataset{k,8}(1+i*binvars:binvars+i*binvars)))};
    end
    mdlsTree{i} = fitctree(cell2mat(dataset(:,1:6)),response);
    mdlsForest{i} = TreeBagger(numTrees,cell2mat(dataset(:,1:6)),response);
    waitbar(i/(N-1),f2,'Training tree')
end
close(f2)
%%
predictorT = microgridTreePredictor(mdlsTree,N,mgrid.ndelta,@mgrid.checkFeasibility);
cSET = microgridSEController(mgrid,N,predictorT);
predictorTB = microgridTreePredictor(mdlsTree,N,mgrid.ndelta,@mgrid.checkFeasibility);
cSETB = microgridSEController(mgrid,N,predictorTB);
%%
[x0,cbuy,csale,cprod,Pload,Pres] = mgrid.representativeSample(N);
input = {x0,cbuy,csale,cprod,Pload,Pres};

tic
[~,r1] = c.control(input);
toc

tic
[~,r2] = cSET.control(input);
toc

tic
[~,r3] = cSETB.control(input);
toc

u1 = r1.V;
u2 = r2.V;
u3 = r3.V;

%%
x1 = zeros(N+1,1);
x1(1) = x0;
x2 = x1;
x3 = x1;

for i = 2:N+1
    [M1,M2,M3] = mgrid.NstepDynamicsMatrices(i-1);
    x1(i) = M1*[u1{1}(1:(i-1)*4);u1{2}(1:(i-1)*4);u1{3}(1:(i-1)*2)] + M2*x0 + M3;
    x2(i) = M1*[u2{1}(1:(i-1)*4);u2{2}(1:(i-1)*4);u2{3}(1:(i-1)*2)] + M2*x0 + M3;
    x3(i) = M1*[u3{1}(1:(i-1)*4);u3{2}(1:(i-1)*4);u3{3}(1:(i-1)*2)] + M2*x0 + M3;
end