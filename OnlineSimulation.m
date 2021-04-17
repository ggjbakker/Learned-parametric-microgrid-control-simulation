addpath('.\classes')
addpath('.\functions\gen5')
%settings
dt = 30;
dtC = 1800;
N = 48;
x0 = 110;

%load data
tv = 0:dt:48*3600-dt;
lt = size(tv,2);
tvh = tv/3600;
dr = 1:48*3600/dt;

cs = zeros(48*3600/dt,1);
cb = zeros(48*3600/dt,1);
cp = zeros(48*3600/dt,1);
cs(1:end/2) = [ones(8*3600/dt,1)*0.18;ones(8*3600/dt,1)*0.5;ones(4*3600/dt,1)*0.45;ones(4*3600/dt,1)*0.32];
cb(1:end/2) = [ones(8*3600/dt,1)*0.55;ones(12*3600/dt,1)*0.775;ones(4*3600/dt,1)*0.55];
cp(1:end/2) = ones(24*3600/dt,1)*0.425;

a = 8000;
b = 14*3600;
Presmax = 700;
Pres = zeros(lt,1);
for i = 1:lt
    ts = tv(i);
    Pres(i) = exp((-(ts-b)^2)/(2*a^2))*Presmax;
end
Pload = load('datasets/Ploaddt15v3.mat');
Pload = Pload.Pload;
Pload = Pload(1:dt/15:end);

%%
%create controller and microgrid classes
mgsim = microgrid(2,dt,0.9,0.9,500,0,-150,150,160,20,-1000,1000,x0);
mgsim3 = microgrid(2,dt,0.9,0.9,500,0,-150,150,160,20,-1000,1000,x0);
mgsim4 = microgrid(2,dt,0.9,0.9,500,0,-150,150,160,20,-1000,1000,x0);
mgctr = microgrid(2,dtC,0.9,0.9,500,0,-150,150,160,20,-1000,1000,x0);

simC1 = simulation(mgsim,dt,cs,cb,cp,Pload,Pres,24*3600/dt,'Regular MPC');
controller = microgridMPCController(mgctr,N);

%simulation
f = waitbar(0,'Simulating');
x = simC1.MLD.x;
for i = 1:lt/2
    cbuyz = cb(i:dtC/dt:i+dtC/dt*(N-1));
    csalez = cs(i:dtC/dt:i+dtC/dt*(N-1));
    cprodz = cp(i:dtC/dt:i+dtC/dt*(N-1));
    Ploadz = Pload(i:dtC/dt:i+dtC/dt*(N-1));
    Presz = Pres(i:dtC/dt:i+dtC/dt*(N-1));
    z = {x,cbuyz,csalez,cprodz,Ploadz,Presz};
    
    tic
    if i == 1
        [u,result] = controller.control(z);
    else
        [u,result] = controller.control(z,result);
    end
    cput = toc;
    
    x = simC1.step(u,cput);
    
    waitbar(i/(lt/2),f,'Simulating')
end
close(f)
simC1.plot()
%%
mgsim3.x = x0;
simC3 = simulation(mgsim3,dt,cs,cb,cp,Pload,Pres,24*3600/dt,'Expression tree parametric MPC');
controller3 = microgridParametricControllerLearned(mgctr,N,predictorTB);

%simulation
f = waitbar(0,'Simulating');
x = simC3.MLD.x;
for i = 1:lt/2
    cbuyz = cb(i:dtC/dt:i+dtC/dt*(N-1));
    csalez = cs(i:dtC/dt:i+dtC/dt*(N-1));
    cprodz = cp(i:dtC/dt:i+dtC/dt*(N-1));
    Ploadz = Pload(i:dtC/dt:i+dtC/dt*(N-1));
    Presz = Pres(i:dtC/dt:i+dtC/dt*(N-1));
    z = {x,cbuyz,csalez,cprodz,Ploadz,Presz};
    
    if i <= dtC/dt
        zm1 = [cbuyz(1);csalez(1);cprodz(1);Ploadz(1);Presz(1)];
        um1 = simC3.u(:,1);
    else
        zm1 = [cb(i-dtC/dt);cs(i-dtC/dt);cp(i-dtC/dt);Pload(i-dtC/dt);Pres(i-dtC/dt)];
        um1 = simC3.u(:,i-dtC/dt);
    end
    
    tic
    if i == 1
        [u,result] = controller3.control(z,zm1,um1);
    else
        [u,result] = controller3.control(z,zm1,um1,result);
    end
    
    cput = toc;  
    x = simC3.step(u,cput);
    
    waitbar(i/(lt/2),f,'Simulating')
end
close(f)
simC3.plot()
%%
mgsim4.x = x0;
simC4 = simulation(mgsim4,dt,cs,cb,cp,Pload,Pres,24*3600/dt,'Hand-designed parametric MPC');
predictor = microgridHeuristicPredictor(N,mgsim4);
controller4 = microgridParametricController(mgctr,N,predictor);

%simulation
f = waitbar(0,'Simulating');
x = simC4.MLD.x;
for i = 1:lt/2
    cbuyz = cb(i:dtC/dt:i+dtC/dt*(N-1));
    csalez = cs(i:dtC/dt:i+dtC/dt*(N-1));
    cprodz = cp(i:dtC/dt:i+dtC/dt*(N-1));
    Ploadz = Pload(i:dtC/dt:i+dtC/dt*(N-1));
    Presz = Pres(i:dtC/dt:i+dtC/dt*(N-1));
    z = {x,cbuyz,csalez,cprodz,Ploadz,Presz};
    
    if i <= dtC/dt
        zm1 = [cbuyz(1);csalez(1);cprodz(1);Ploadz(1);Presz(1)];
        um1 = simC4.u(:,1);
    else
        zm1 = [cb(i-dtC/dt);cs(i-dtC/dt);cp(i-dtC/dt);Pload(i-dtC/dt);Pres(i-dtC/dt)];
        um1 = simC4.u(:,i-dtC/dt);
    end
    
    tic
    if i == 1
        [u,result] = controller4.control(z,zm1,um1);
    else
        [u,result] = controller4.control(z,zm1,um1,result);
    end
    
    cput = toc;  
    x = simC4.step(u,cput);
    
    waitbar(i/(lt/2),f,'Simulating')
end
close(f)
simC4.plot()

run plotfigureonline.m