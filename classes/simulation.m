classdef simulation < handle
    properties
        MLD
        dt
        csale
        cbuy
        cprod
        Pload
        Pres
        Tmax
        totalcost
        importcost
        exportcost
        gencost
        computationtime
        k
        kmax
        u
        x
        controllername
    end
    methods
        function obj = simulation(MLD,dt,cs,cb,cp,Pload,Pres,kmax,controllername)
            obj.MLD = MLD;
            obj.Pload = Pload;
            obj.Pres = Pres;
            obj.dt = dt;
            obj.csale = cs;
            obj.cbuy = cb;
            obj.cprod = cp;
            obj.kmax = kmax;
            obj.controllername = controllername;
            
            
            obj.k = 0;
            obj.totalcost = zeros(1,kmax);
            obj.importcost = zeros(1,kmax);
            obj.exportcost = zeros(1,kmax);
            obj.gencost = zeros(1,kmax);
            obj.computationtime = zeros(1,kmax);
            obj.u = zeros(MLD.nu,kmax);
            obj.x = zeros(MLD.nx,kmax);
        end
        function xt = step(obj,u,walltime)
            if sum(u(2:obj.MLD.nu)) - u(1) + obj.Pres(obj.k+1) - obj.Pload(obj.k+1) > 0.001
                error('Power balance violation')
            end
            xt = obj.MLD.step(u);
            obj.k = obj.k+1;
            obj.x(obj.k) = xt;
            obj.u(:,obj.k) = u;
            obj.computationtime(obj.k) = walltime;
            if u(4) > 0
                obj.importcost(obj.k) = u(4)*obj.cbuy(obj.k);
                obj.exportcost(obj.k) = 0;
                obj.gencost(obj.k) = sum(u(2:end-1))*obj.cprod(obj.k);
                obj.totalcost(obj.k) = sum(u(2:end-1))*obj.cprod(obj.k) + u(4)*obj.cbuy(obj.k);
            else
                obj.importcost(obj.k) = 0;
                obj.exportcost(obj.k) = u(4)*obj.csale(obj.k);
                obj.gencost(obj.k) = sum(u(2:end-1))*obj.cprod(obj.k);
                obj.totalcost(obj.k) = sum(u(2:end-1))*obj.cprod(obj.k) + u(4)*obj.csale(obj.k);
            end
        end
        function plot(obj)
            dr = 1:obj.k;
            tv = 0:obj.dt:obj.kmax*obj.dt;
            tvh = tv/3600;
            
            Pgrid = obj.u(end,:);
            Pb = obj.u(1,:);
            Ptot = sum(obj.u(2:end-1,:),1);
            Pgen1 = obj.u(2,:);
            Pgen2 = obj.u(3,:);
    
            fg = figure;
            fg.Position =  [100 100 800 400];
            plot(tvh(dr),obj.Pres(dr),tvh(dr),obj.Pload(dr),tvh(dr),Pgrid(dr),tvh(dr),Pb(dr),tvh(dr),Ptot(dr))
            title(obj.controllername,'Interpreter','latex','FontSize',16)
            ylabel({'$P(\mathrm{kW})$'},'Interpreter','latex','FontSize',14)
            xlabel({'t (hours)'},'Interpreter','latex','FontSize',16)
            legend({'$P_\mathrm{res}$','$P_\mathrm{load}$','$P_\mathrm{grid}$','$P_\mathrm{ess}$','$P_\mathrm{gen,tot}$'},'Interpreter','latex','FontSize',12,'Location','SouthEast')
            fg = figure;
            fg.Position =  [100 100 800 400];
            plot(tvh(dr),obj.x(dr))
            title(obj.controllername,'Interpreter','latex','FontSize',16)
            xlabel({'t (hours)'},'Interpreter','latex','FontSize',16)
            ylabel({'Energy ($\mathrm{kWh}$)'},'Interpreter','latex','FontSize',14)
            legend({'$x_\mathrm{ess}$'},'Interpreter','latex','FontSize',12,'location','northwest')
        end
    end
end