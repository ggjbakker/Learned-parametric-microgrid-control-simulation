dr = 1:24*3600/dt;
fg = figure;
fg.Position =  [100 100 800 400];
plot(tvh(dr),Pres(dr),tvh(dr),Pload(dr))
title({'Case study scenario'},'Interpreter','latex','FontSize',16)
ylabel({'$P(\mathrm{kW})$'},'Interpreter','latex','FontSize',14)
xlabel({'t (hours)'},'Interpreter','latex','FontSize',16)
legend({'$P_\mathrm{res}$','$P_\mathrm{load}$'},'Interpreter','latex','FontSize',16)
%%
fg = figure;
fg.Position =  [100 100 800 400];
plot(tvh(dr),cp(dr),tvh(dr),cb(dr),tvh(dr),cs(dr))
title({'Case study scenario'},'Interpreter','latex','FontSize',16)
ylabel({'cost'},'Interpreter','latex','FontSize',14)
xlabel({'t (hours)'},'Interpreter','latex','FontSize',16)
legend({'$c_\mathrm{gen}$','$c_\mathrm{import}$','$c_\mathrm{export}$'},'Interpreter','latex','FontSize',16)
ylim([0 1])
%%
fg = figure;
fg.Position =  [100 100 600 400];
p = uipanel('Parent',fg,'BorderType','none'); 
p.Title = 'Computation time'; 
p.TitlePosition = 'centertop'; 
p.FontSize = 12;
p.BackgroundColor = [1 1 1];
sp = subplot(1,3,1,'Parent',p);
boxplot([simC1.computationtime'],'Labels',{'Regular MPC'})
ylabel('time (s)')
s= subplot(1,3,2,'Parent',p);
boxplot([simC3.computationtime'],'Labels',{'Expression tree MPC'})
s= subplot(1,3,3,'Parent',p);
boxplot([simC4.computationtime'],'Labels',{'Hand designed MPC'})