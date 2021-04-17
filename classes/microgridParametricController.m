classdef microgridParametricController < MLDController
    properties
        N %horizon of MPC controller
        MLD %MLD model
        
        %matrices for updating cost
        mz
        mu1
        mu2
        
        %matrices for updating bounds
        bupdate1
        bupdate2
        bupdate3
        
        M1
        M2
        
        %optimization matrices
        A1
        A2
        A3
        A4
        
        b1
        b2
        b3
        b4
        b5
        
        c

        %predictor
        predictor
        
        %parametric variables      
        prob
        
        options
    end
    methods
        %%
        %Constructor function
        %
        function obj = microgridParametricController(MLD,N,predictor)
            
            obj.N = N;
            obj.MLD = MLD;
            
            obj.setDynamicMatrices();
            obj.setCostMatrices();
            obj.createProb();
            
            obj.predictor = predictor;
            
            obj.options = struct;
            obj.prob.LargeScale = 1;
%             obj.prob.NumDiff = 6;
%             obj.prob.ConsDiff = 6;
            
        end     
        %%
        %
        %
        function [u,result] = control(obj,z,zm1,u0,theta0) 
            delta = obj.predictor.predict(z);
            obj.updateBounds(z);
            obj.updateCost(z);
            
            obj.prob.user.u0 = u0;
            obj.prob.user.x0 = z{1};
            obj.prob.user.delta = delta;

            obj.prob.user.cprod = z{4};
            obj.prob.user.csale = z{3};
            obj.prob.user.Pload = z{5};
            obj.prob.user.Pres = z{6};
            obj.prob.user.cbuy = z{2};
            obj.prob.user.Pload0 = zm1(4);
            obj.prob.user.Pres0 =zm1(5);

            obj.prob.user.cs = obj.c;
            
            obj.prob = modify_c_U(obj.prob,10^10*ones(obj.MLD.nx*obj.N*2+obj.MLD.nconstraints*obj.N+obj.N+obj.N*(obj.MLD.nu+obj.MLD.nz),1));
            obj.prob = modify_c_L(obj.prob,[-1e4*ones(size(obj.b1));-1e4*ones(size(obj.b2));obj.b3;obj.b4]);
            obj.prob = modify_c_U(obj.prob,[obj.b1;obj.b2;obj.b3;obj.b5]);
            if nargin > 4
                obj.prob = WarmDefSOL('snopt', obj.prob, theta0);
            end

            result = tomRun('snopt',obj.prob,0);
            V = conPar_V(result.x_k,obj.prob);
            result.V = V;
            
            u = V{1}(1:4);
            if abs(u(2)) < 1
                u(2) = 0;
            end
            if abs(u(3)) < 1
                u(3) = 0;
            end
            if u(2) < obj.MLD.Pdis_l && u(2) >= 1
                u(2) = obj.MLD.Pdis_l;
            end
            if u(3) < obj.MLD.Pdis_l && u(3) >= 1
                u(3) = obj.MLD.Pdis_l;
            end
            
            u = min(obj.MLD.uub,max(obj.MLD.ulb,u));
            if z{1} <  obj.MLD.xub - ((obj.MLD.xub-obj.MLD.xlb)*0.95 + obj.MLD.xlb) && u(1) < 0
                u(1) = 0;
            end
            if z{1} > obj.MLD.xub*0.95 && u(1) > 0
                u(1) = 0;
            end

            u(4) = u(4) - z{6}(1) - sum(u(2:4)) + z{5}(1) + u(1);
        end
        %%
        %
        %
        function setDynamicMatrices(obj)
            [P1,P2,P3] = obj.MLD.statePredictionMatrices(obj.N);
            [F1,F2,F3] = obj.MLD.NstepConstraintMatrices(obj.N);
            obj.bupdate1 = F2;
            obj.M1 = F3;
            obj.M2 = P2;
            
            %           
            obj.A1 = F1;
            obj.b1 = zeros(obj.MLD.nconstraints*obj.N,1);
            
            obj.A2 = [P1;-P1];
            obj.b2 = zeros(obj.MLD.nx*obj.N*2,1);
            
            m = -1*ones(1,obj.MLD.nu);
            m(1) = 1;
            obj.A3 = [kron(eye(obj.N),m),zeros(obj.N,obj.N*(obj.MLD.ndelta+obj.MLD.nz))];
            obj.b3 = zeros(obj.N,1);
            
            m = [eye(obj.MLD.nu*obj.N),zeros(obj.MLD.nu*obj.N,obj.MLD.ndelta*obj.N),zeros(obj.MLD.nu*obj.N,obj.MLD.nz*obj.N) ;
                    zeros(obj.MLD.nz*obj.N,obj.MLD.nu*obj.N),zeros(obj.MLD.nz*obj.N,obj.MLD.ndelta*obj.N),eye(obj.MLD.nz*obj.N)];
            obj.A4 = m;
            obj.b4 = [kron(ones(obj.N,1),obj.MLD.ulb);kron(ones(obj.N,1),obj.MLD.zlb)];
            obj.b5 = [kron(ones(obj.N,1),obj.MLD.uub);kron(ones(obj.N,1),obj.MLD.zub)];
            
            %
            obj.bupdate2 = ones(obj.MLD.nx*obj.N,1)*obj.MLD.xub - P3;
            obj.bupdate3 = -ones(obj.MLD.nx*obj.N,1)*obj.MLD.xlb + P3;
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
        %
        %
        %
        function createProb(obj)
            obj.prob = conAssign('conPar_f','conPar_g',[],[],[],[],'parMPC',zeros(10,1),[],'-1E10',[],[],[],'conPar_c',[],[],[],-10^10*ones(obj.MLD.nx*obj.N*2+obj.MLD.nconstraints*obj.N+obj.N+obj.N*(obj.MLD.nu+obj.MLD.nz),1),10^10*ones(obj.MLD.nx*obj.N*2+obj.MLD.nconstraints*obj.N+obj.N+obj.N*(obj.MLD.nu+obj.MLD.nz),1));

            obj.prob.user.mA = obj.MLD.A;
            obj.prob.user.mB1 = obj.MLD.B1;
            obj.prob.user.mB2 = obj.MLD.B2;
            obj.prob.user.mB3 = obj.MLD.B3;

            obj.prob.user.xub = obj.MLD.xub;
            obj.prob.user.xlb = obj.MLD.xlb;

             obj.prob.user.A1 = obj.A1;
             obj.prob.user.A2 = obj.A2;
             obj.prob.user.A3 = obj.A3;
             obj.prob.user.A4 = obj.A4;
             obj.prob.user.N = obj.N;
        end
        %%
        %Calculates new bounds
        %
        function updateBounds(obj,z)
            %parameters
            x0 = z{1};
            Pload = z{5};
            Pres = z{6};
            
            %update bounds
            obj.b1 = obj.bupdate1+obj.M1*x0;
            
            obj.b2 = [obj.bupdate2 - obj.M2*x0;obj.bupdate3 + obj.M2*x0];
            
            obj.b3 = Pres-Pload;
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