classdef microgridMPCControllerTOMLAB < MLDController
    properties
        N %horizon of QMPC controller
        microgrid
        prob
        
        %optimization matrices
        mz
        mu1
        mu2
        c
        
        A
        b_L
        b_U
        
        M2
        M3
        F2
        F3
        
        x_L
        x_U
        
        intVars
    end
    methods
        %%
        %Constructor function
        %
        function obj = microgridMPCControllerTOMLAB(microgrid,N)     
            obj.N = N;
            obj.microgrid = microgrid;
            
            obj.setMatrices();
            obj.prob = mipAssign(obj.c,obj.A,obj.b_L,obj.b_U,obj.x_L,obj.x_U,[],'microgridMILPMPC',[],[],obj.intVars);
            obj.prob.Solver.Alg = 0;
            obj.prob.optParam.MaxIter = 50000;
            obj.prob.Threads = 6;
        end     
        %%
        %
        %
        function [u,result] = control(obj,z)
            obj.updateBounds(z);
            obj.prob = modify_b_L(obj.prob,obj.b_L);
            obj.prob = modify_b_U(obj.prob,obj.b_U);
            obj.prob = modify_c(obj.prob,obj.c);
            result = tomRun('milpsolve', obj.prob, 0);
            u = {};
            u{1} = result.x_k(1:obj.microgrid.nu*obj.N);
            u{2} = result.x_k(obj.microgrid.nu*obj.N+1:(obj.microgrid.nu+obj.microgrid.ndelta)*obj.N);
            u{3} = result.x_k((obj.microgrid.nu+obj.microgrid.ndelta)*obj.N+1:end);
            result.V = u;
            u = u{1}(1:4);
        end
        %%
        %MILP controller
        %
        function setMatrices(obj)
            N = obj.N;
            %cost matrices
            obj.mz = zeros(2*N,N);
            obj.mz(2:2:end,:) = eye(N);
            obj.mu1 = zeros(obj.microgrid.nu*N,N);
            for i = 0:N-1
                for j = 0:obj.microgrid.ndis-1
                    obj.mu1(2+obj.microgrid.nu*i+j,i+1) = 1;
                end
            end
            obj.mu2 = zeros(obj.microgrid.nu*N,N);
            obj.mu2(2+obj.microgrid.ndis:obj.microgrid.nu:end,:) = eye(N);
            obj.c = zeros((obj.microgrid.nu+obj.microgrid.ndelta+obj.microgrid.nz)*N,1);
            
            %MLDdynamic constraints
            M1 = zeros(obj.microgrid.nx*N,(obj.microgrid.nu+obj.microgrid.ndelta+obj.microgrid.nz)*N);
            obj.M2 = zeros(obj.microgrid.nx*N,obj.microgrid.nx);
            obj.M3 = zeros(obj.microgrid.nx*N,1);
            
            nx = obj.microgrid.nx;
            di = obj.microgrid.nu*N;
            dz = (obj.microgrid.nu+obj.microgrid.ndelta)*N;
            for i = 1:N
                if i == 1
                    M1(1:nx,1:obj.microgrid.nu) = obj.microgrid.B1;
                    
                    M1(1:nx,di+1:di+obj.microgrid.ndelta) = obj.microgrid.B2;
                    
                    M1(1:nx,dz+1:dz+obj.microgrid.nz) = obj.microgrid.B3;
                    
                    obj.M3(1:nx) = obj.microgrid.B4;
                else
                    M1(1+(i-1)*nx:i*nx,1+obj.microgrid.nu:obj.microgrid.nu*N) = M1(1+(i-2)*nx:(i-1)*nx,1:obj.microgrid.nu*(N-1));
                    M1(1+(i-1)*nx:i*nx,1:obj.microgrid.nu) = obj.microgrid.A*M1(1+(i-1)*nx:i*nx,1+obj.microgrid.nu:2*obj.microgrid.nu);

                    M1(1+(i-1)*nx:i*nx,1+di+obj.microgrid.ndelta:di+obj.microgrid.ndelta*N) = M1(1+(i-2)*nx:(i-1)*nx,1+di:di+obj.microgrid.ndelta*(N-1));
                    M1(1+(i-1)*nx:i*nx,1+di:di+obj.microgrid.ndelta) = obj.microgrid.A*M1(1+(i-1)*nx:i*nx,1+di+obj.microgrid.ndelta:di+2*obj.microgrid.ndelta);

                    M1(1+(i-1)*nx:i*nx,1+dz+obj.microgrid.nz:dz+obj.microgrid.nz*N) = M1(1+(i-2)*nx:(i-1)*nx,1+dz:dz+obj.microgrid.nz*(N-1));
                    M1(1+(i-1)*nx:i*nx,1+dz:dz+obj.microgrid.nz) = obj.microgrid.A*M1(1+(i-1)*nx:i*nx,1+dz+obj.microgrid.nz:dz+2*obj.microgrid.nz);

                    obj.M3(1+(i-1)*nx:i*nx) = obj.M3(1+(i-2)*nx:(i-1)*nx,:) + obj.microgrid.A^i*obj.microgrid.B4;
                end
                
                obj.M2(1+(i-1)*nx:i*nx,:) = obj.microgrid.A^i;
            end
            
            [F1,F2,F3,~,~,~] = obj.microgrid.NstepMatrices(i);
            obj.F2 = F2;
            obj.F3 = F3;
            
            %
            noc = obj.microgrid.nconstraints*N+obj.microgrid.nx*N+N;
            obj.A = zeros(noc,(obj.microgrid.nu+obj.microgrid.ndelta+obj.microgrid.nz)*N);
            obj.b_U = zeros(noc,1);
            obj.b_L =  zeros(noc,1);
            
            cend = obj.microgrid.nconstraints*N;
            obj.A(1:cend,:) = F1;
            cstart = cend;
            
            cend = cstart+obj.microgrid.nx*N;
            obj.A(cstart+1:cend,:) = M1;
            cstart = cend;
            
            cend = cstart+N;
            m = -1*ones(1,obj.microgrid.nu);
            m(1) = 1;
            obj.A(cstart+1:cend,1:obj.microgrid.nu*N) = kron(eye(N),m);
            
            %
            obj.x_L = [kron(ones(N,1),obj.microgrid.ulb);zeros(N*obj.microgrid.ndelta,1);kron(ones(N,1),obj.microgrid.zlb)];
            obj.x_U = [kron(ones(N,1),obj.microgrid.uub);ones(N*obj.microgrid.ndelta,1);kron(ones(N,1),obj.microgrid.zub)];
            obj.intVars = [zeros(N*obj.microgrid.nu,1);ones(N*obj.microgrid.ndelta,1);zeros(N*obj.microgrid.nz,1)];
            
            %
            obj.A = sparse(obj.A);
        end
        %%
        %Calculates new bounds and cost function
        %
        function updateBounds(obj,z)
            
            %parameters
            x0 = z{1};
            cbuy = z{2};
            csale = z{3};
            cprod = z{4};
            Pload = z{5};
            Pres = z{6};
            
            %cost matrix
            cu = (obj.mu1*cprod)+(obj.mu2*cbuy);
            cd = zeros(obj.microgrid.ndelta*obj.N,1);
            cz = obj.mz*(csale-cbuy);
            obj.c = [cu;cd;cz]; 
            
            %A matrix
            cstart = 0;
            cend = obj.microgrid.nconstraints*obj.N;
            obj.b_U(cstart+1:cend) = obj.F2+obj.F3*x0;
            obj.b_L(cstart+1:cend) = -inf;
            cstart = cend;
            
            cend = cstart+obj.microgrid.nx*obj.N;
            obj.b_L(cstart+1:cend) = ones(obj.microgrid.nx*obj.N,1)*obj.microgrid.xlb - obj.M3 - obj.M2*x0+obj.microgrid.epsilon;
            obj.b_U(cstart+1:cend) = ones(obj.microgrid.nx*obj.N,1)*obj.microgrid.xub + obj.M3 + obj.M2*x0;
            cstart = cend;
            
            cend = cstart+obj.N;
            obj.b_L(cstart+1:cend) = Pres-Pload;
            obj.b_U(cstart+1:cend) = Pres-Pload;
        end
    end
end