classdef microgrid < MLDModel
    properties
        ndis %number of distributed generators
        Pdis_l %lower limit of power production for distributed generators
        sT %sampling time
        rprfilter %filter for creating random samples
    end
    methods
        %Constructor method
        %
        %
        %
        %
        function obj = microgrid(ndis,sT,eta_c,eta_d,x_h,x_l,Pb_l,Pb_h,Pdis_h,Pdis_l,Pgrid_l,Pgrid_h,x0)
            if nargin < 3
                %Default values for microgrid parameters
                eta_c = 0.92; %charging efficiency [-]
                eta_d = 0.87; %discharging efficiency [-]
                x_l = 0;
                x_h = 36; %upper energy storage limit [kWh]
                x0 = 12;
                Pb_l = -3; %max discharge rate [kW]
                Pb_h = 2; %max scharge rate [kW]

                Pdis_h = 6;
                Pdis_l = 4;
                Pgrid_l = -15;
                Pgrid_h = 15;
            end
            
            epsilon = 10^-6;
            %bounds
            xlb = [x_l];
            xub = [x_h];
            zlb = [Pb_l];  
            zub = [Pb_h];
            ulb = [Pb_l];
            uub = [Pb_h];
            
            %state names
            xname = ["xb"];
            uname = ["Pb"];
            zname = ["db*Pb"];
            dname = ["db"];
            
            obj = obj@MLDModel([],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],0,"x0",x0);
            obj.addInputs(1,ulb,uub,uname);
            obj.addDeltas(1,dname);
            obj.addStates(1,xlb,xub,[],[],[],[],[],[],[],xname)
            obj.addBinaryEquivalence(0,1,0,0,dname)
            obj.addBinaryContinuousMultiplication([],1,[],[],dname,"uname",uname,"variablename",zname)
            obj.A = 1;
            obj.B1 = sT/3600*eta_c;
            obj.B3 = sT/3600*(1/eta_d - eta_c);
            
            %Add generators
            unames = {};
            dnames = {};
            ub = [];
            lb = [];
            for i = 1:ndis
                unames = [unames,"Pdis"+string(i)];
                dnames = [dnames,"d"+string(i)+"on"];
                ub = [ub;Pdis_h];
                lb = [lb;0];
            end
            obj.addInputs(ndis,lb,ub,unames);
            obj.addDeltas(ndis,dnames);
            for i = 1:ndis
                obj.addConstraints(0,-1,Pdis_l,0,0,"uname",{"Pdis"+string(i)},"dname",{"d"+string(i)+"on"});
                obj.addConstraints(0,1,-Pdis_h,0,0,"uname",{"Pdis"+string(i)},"dname",{"d"+string(i)+"on"});
            end
            
            %Add grid
            obj.addInputs(1,Pgrid_l,Pgrid_h,{"Pgrid"});
            obj.addDeltas(1,{"dgrid"});
            obj.addBinaryEquivalence(0,1,0,0,"dgrid","uname","Pgrid");
            
            obj.addBinaryContinuousMultiplication(0,1,0,0,"dgrid","uname","Pgrid","variablename","dgrid*Pgrid");
            
            obj.ndis = ndis;
            obj.Pdis_l = Pdis_l;
            
            
            %filter used to create representative samples
            obj.rprfilter = designfilt('lowpassiir','FilterOrder',3, ...
             'PassbandFrequency',2.2e-5,'PassbandRipple',0.2, ...
             'SampleRate',1/sT);
            obj.sT = sT;
        end
        %%
        %
        %
        %
        function [feasibility, corrected] = checkFeasibility(obj,deltas,Pload,Pres)
            feasibility = 1;
            corrected = deltas;
            if (sum(deltas(2:1+obj.ndis)*obj.uub(2))-deltas(1)*obj.ulb(1) < Pload-Pres) && (deltas(end) == 1)
                feasibility = 0;
                corrected(end) = 0;
            end

            if ((Pres-Pload)+ sum(deltas(2:1+obj.ndis)*obj.Pdis_l) > (1-deltas(1))*obj.uub(1)) && (deltas(end) == 0)
                feasibility = 0;
                corrected(end) = 1;
            end

             if 0 < (Pload - Pres) && (Pload-Pres) < sum(deltas(2:1+obj.ndis)*obj.Pdis_l)-(1-deltas(1))*obj.uub(1) && (deltas(end) == 0)
                 feasibility = 0;
                 corrected(end) = 1;
             end
        end
        %%
        %
        %
        %
        function Zn = representativeSet(obj,Np,N)
            Zn = cell(N,6);
            for i = 1:N
                [x0,cbuy,csale,cprod,Pload,Pres] = obj.representativeSample(Np);
                Zn(i,:) = {x0,cbuy,csale,cprod,Pload,Pres};
            end
        end
        %%
        %
        %
        %
        function [x0,cbuy,csale,cprod,Pload,Pres] = representativeSample(obj,Np)
            x0 = rand(1)*(obj.xub-obj.xlb)+obj.xlb;
            
            cbuy = filter(obj.rprfilter,rand(Np,1)-0.5)*60/obj.sT^0.5+1;
            cprod = filter(obj.rprfilter,rand(Np,1)-0.5)*60/obj.sT^0.5+0.75;
            csale = cbuy*0.8;
            Pload = filter(obj.rprfilter,rand(Np,1)-0.5)*500/obj.sT^0.5+15;
            Pload = Pload*40;
            Pres = filter(obj.rprfilter,rand(Np,1)-0.5)*1500/obj.sT^0.5+10;
            Pres = Pres*30;
            
        end
        %%
        %
        %
        function plot(obj,N,X,cbuy,csale,cprod,Pload,Pres,V)
            figure
            plot(0:N-1,cbuy,0:N-1,cprod,0:N-1,csale)
            legend('cbuy','cprod','csale')
            title('Energy prices')
            axis([0 N 0 2])
            
            figure
            hold on
            axis([0 N obj.ulb(4) obj.uub(4)])
            plot(0:N-1,V(1:2+obj.ndis:(N-1)*(2+obj.ndis)+1))
            Legend = cell(4+obj.ndis,1);
            Legend{1} = 'Pbattery';
            for i = 1:obj.ndis
                plot(0:N-1,V(1+i:2+obj.ndis:(N-1)*(2+obj.ndis)+1+i))
                Legend{i+1}=strcat('Pdis', num2str(i));
            end
            plot(0:N-1,V(1+obj.ndis+1:2+obj.ndis:(N-1)*(2+obj.ndis)+2+obj.ndis))
            Legend{obj.ndis+2} = 'Pgrid';
            plot(0:N-1,Pload)
            Legend{obj.ndis+3} = 'Pload';
            plot(0:N-1,Pres)
            Legend{obj.ndis+4} = 'Pres';
            
            legend(Legend)
            title('Energy production')
            
            figure
            hold on
            axis([0 N obj.xlb(1) obj.xub(1)])
            plot(0:N,X)
            title('Battery charge')
            legend('xb')
        end
        %%
        %test function to check optimization
        %
        function test(obj,N,opts)
            if nargin < 3
                opts = obj.defopts;
            end
            %inputs and states
            xest = sdpvar((N+1)*obj.nx,1);
            
            U = sdpvar(N*obj.nu,1); 
            DELTA = binvar(N*obj.ndelta,1);
            Z = sdpvar(N*obj.nz,1);
            
            V = [U;DELTA;Z];
            
            %controller parameters
            [x0,cbuy,csale,cprod,Pload,Pres] = obj.representativeSample(N);
            constraints = [xest(1) == x0];
            for i = 2:N+1
                [F1,F2,F3,M1,M2,M3] = obj.NstepMatrices(i-1);
                constraints = [constraints,obj.xlb <= M1*[U(1:(i-1)*4);DELTA(1:(i-1)*4);Z(1:(i-1)*2)] + M2 * x0 + M3 <= obj.xub]; %#ok<*AGROW>
                constraints = [constraints, xest(1+(i-1)*obj.nx:i*obj.nx) == (M1*[U(1:(i-1)*4);DELTA(1:(i-1)*4);Z(1:(i-1)*2)] + M2 * x0 + M3)]; %could remove
            end
            
            constraints = [constraints,F1*V <= F2 + F3*x0];
            constraints = [constraints, kron(ones(N,1),obj.ulb) <= U <= kron(ones(N,1),obj.uub)]; %#ok<*CHAIN>
            constraints = [constraints, kron(ones(N,1),obj.zlb) <= Z <= kron(ones(N,1),obj.zub)];      
            m = -1*ones(1,obj.nu);
            m(1) = 1;
            constraints = [constraints, kron(eye(N),m)*U == Pres-Pload];
            
            mz = zeros(2*N,N);
            mz(2:2:end,:) = eye(N);
            Cz = (mz*(csale-cbuy))';
            
            mu1 = zeros(obj.nu*N,N);
            for i = 0:N-1
                for j = 0:obj.ndis-1
                    mu1(2+obj.nu*i+j,i+1) = 1;
                end
            end
            mu2 = zeros(obj.nu*N,N);
            mu2(2+obj.ndis:obj.nu:end,:) = eye(N);

            Cu = (mu1*cprod)'+(mu2*cbuy)';
            
            Cd = ones(1,obj.ndelta*N);
            C = [Cu,Cd,Cz];
            objective = C*V;
            
            result = optimize(constraints,objective,opts);
            
            obj.plot(N,value(xest),cbuy,csale,cprod,Pload,Pres,value(V));
        end
    end
end