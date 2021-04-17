classdef microgridMPCController < MLDController
    properties
        N %horizon of MPC controller
        MLD %MLD model
        prob
        params
        
        %matrices for updating cost
        mz
        mu1
        mu2
        
        %matrices for updating bounds
        b1
        b2
        b3
        M1
        M2
        
        %optimization matrices
        A
        b
        c
        x_lb
        intVars
        sense
    end
    methods
        %%
        %Constructor function
        %
        function obj = microgridMPCController(MLD,N)     
            obj.N = N;
            obj.MLD = MLD;
            
            obj.setDynamicMatrices();
            obj.setCostMatrices();
            
            obj.prob = {};
            obj.prob.A = sparse(obj.A);
            obj.prob.obj = obj.c';
            obj.prob.rhs = obj.b;
            obj.prob.lb = obj.x_lb;
            obj.prob.sense = obj.sense;
            obj.prob.vtype = char(obj.intVars*'B'-1*(obj.intVars-1)*'C');
            
            obj.params.outputflag = 0;
        end     
        %%
        %
        %
        function [u,result] = control(obj,z,result)
            obj.updateBounds(z);
            obj.updateCost(z);
            obj.prob.obj = obj.c';
            obj.prob.rhs = obj.b;
            if nargin > 2
                obj.prob.start = result.x;
            end
            result = gurobi(obj.prob, obj.params);
            
            u = {};
            u{1} = result.x(1:obj.MLD.nu*obj.N);
            u{2} = result.x(obj.MLD.nu*obj.N+1:(obj.MLD.nu+obj.MLD.ndelta)*obj.N);
            u{3} = result.x((obj.MLD.nu+obj.MLD.ndelta)*obj.N+1:end);
            result.V = u;
            u = u{1}(1:obj.MLD.nu);
        end
        %%
        %
        %
        function setDynamicMatrices(obj)
            [P1,P2,P3] = obj.MLD.statePredictionMatrices(obj.N);
            [F1,F2,F3] = obj.MLD.NstepConstraintMatrices(obj.N);
            obj.b1 = F2;
            obj.M1 = F3;
            obj.M2 = P2;
            
            %
            noc = obj.MLD.nconstraints*obj.N+obj.MLD.nx*obj.N*2+obj.N+(obj.MLD.nu+obj.MLD.nz)*obj.N*2;
            obj.A = zeros(noc,(obj.MLD.nu+obj.MLD.ndelta+obj.MLD.nz)*obj.N);
            obj.b = zeros(noc,1);
            obj.sense = char(ones(noc,1)*'<');
            
            noc = obj.MLD.nconstraints*obj.N;
            obj.A(1:noc,:) = F1;
            cstart = noc;
            
            noc = obj.MLD.nx*obj.N;
            obj.A(cstart+1:cstart+noc*2,:) = [P1;-P1];
            obj.sense(cstart+1:cstart+noc) = char(ones(noc,1)*'<');
            cstart = cstart+noc*2;
            
            noc = obj.N;
            m = -1*ones(1,obj.MLD.nu);
            m(1) = 1;
            obj.A(cstart+1:cstart+noc,1:obj.MLD.nu*obj.N) = kron(eye(obj.N),m);
            obj.sense(cstart+1:cstart+noc) = char(ones(noc,1)*'=');
            cstart = cstart+noc;
            
            noc = (obj.MLD.nu+obj.MLD.nz)*obj.N;
            m = [eye(obj.MLD.nu*obj.N),zeros(obj.MLD.nu*obj.N,obj.MLD.ndelta*obj.N),zeros(obj.MLD.nu*obj.N,obj.MLD.nz*obj.N) ;
                    zeros(obj.MLD.nz*obj.N,obj.MLD.nu*obj.N),zeros(obj.MLD.nz*obj.N,obj.MLD.ndelta*obj.N),eye(obj.MLD.nz*obj.N)];
            obj.A(cstart+1:cstart+2*noc,:) = [m;m];
            obj.sense(cstart+1:cstart+noc) = char(ones(noc,1)*'>');
            obj.b(cstart+1:cstart+noc) = [kron(ones(obj.N,1),obj.MLD.ulb);kron(ones(obj.N,1),obj.MLD.zlb)];
            obj.b(cstart+noc+1:cstart+noc*2) = [kron(ones(obj.N,1),obj.MLD.uub);kron(ones(obj.N,1),obj.MLD.zub)];
            
            %
            obj.x_lb = [kron(ones(obj.N,1),obj.MLD.ulb);zeros(obj.N*obj.MLD.ndelta,1);kron(ones(obj.N,1),obj.MLD.zlb)];
            obj.intVars = [zeros(obj.N*obj.MLD.nu,1);ones(obj.N*obj.MLD.ndelta,1);zeros(obj.N*obj.MLD.nz,1)];
            
            %
            obj.b2 = ones(obj.MLD.nx*obj.N,1)*obj.MLD.xub - P3;
            obj.b3 = -ones(obj.MLD.nx*obj.N,1)*obj.MLD.xlb + P3;
        end
        %%
        %
        %
        function setCostMatrices(obj)
            %cost matrices
            obj.mz = zeros(2*obj.N,obj.N);
            obj.mz(2:2:end,:) = eye(obj.N);
            obj.mu1 = zeros(obj.MLD.nu*obj.N,obj.N);
            for i = 0:obj.N-1
                for j = 0:obj.MLD.ndis-1
                    obj.mu1(2+obj.MLD.nu*i+j,i+1) = 1;
                end
            end
            obj.mu2 = zeros(obj.MLD.nu*obj.N,obj.N);
            obj.mu2(2+obj.MLD.ndis:obj.MLD.nu:end,:) = eye(obj.N);
            obj.c = zeros((obj.MLD.nu+obj.MLD.ndelta+obj.MLD.nz)*obj.N,1);
        end
        %%
        %Calculates new bounds and cost function
        %
        function updateBounds(obj,z)      
            %parameters
            x0 = z{1};
            Pload = z{5};
            Pres = z{6};
            
            %update bounds
            cstart = 0;
            noc = obj.MLD.nconstraints*obj.N;
            obj.b(cstart+1:noc) = obj.b1+obj.M1*x0;
            cstart = noc;
            
            noc = obj.MLD.nx*obj.N;
            obj.b(cstart+1:cstart+noc) = obj.b2 - obj.M2*x0;
            obj.b(cstart+noc+1:cstart+noc*2) = obj.b3 + obj.M2*x0;
            cstart = cstart+2*noc;
            
            noc = obj.N;
            obj.b(cstart+1:cstart+noc) = Pres-Pload;
        end
        %%
        %
        %
        function updateCost(obj,z)
            %parameters
            cbuy = z{2};
            csale = z{3};
            cprod = z{4};
            %update cost
            cu = (obj.mu1*cprod)+(obj.mu2*cbuy);
            cd = zeros(obj.MLD.ndelta*obj.N,1);
            cz = obj.mz*(csale-cbuy);
            obj.c = [cu;cd;cz]; 
        end
    end
end