close all
addpath('./classes')
load('C:\Users\geert\OneDrive\Systems and control\Afstudeer simulatie\datasets\Genpar5\checkpoint.mat')
kp = 65;
gvec = 1:150:kp*150;
figure
p = plot(gvec,best2(1:kp),'color',[0.2 0.2 0.2]);
p(1).LineWidth=1.2;
title({'Fitness best individual $l_{\mathrm{gen}}$'},'Interpreter','latex','FontSize',16)
ylabel({'$\mathrm{f}_i$'},'Interpreter','latex','FontSize',16)
xlabel({'Generation'},'Interpreter','latex','FontSize',16)
figure
p = plot(gvec,best4(1:kp),'color',[0.2 0.2 0.2]);
p(1).LineWidth=1.2;
title({'Fitness best individual $l_{\mathrm{grid}}$'},'Interpreter','latex','FontSize',16)
ylabel({'$\mathrm{f}_i$'},'Interpreter','latex','FontSize',16)
xlabel({'Generation'},'Interpreter','latex','FontSize',16)