classdef MLDModel < handle
    properties
        %Dynamic matrices
        %   x(k+1) = A*x + B1*u + B2*delta + B3*z + B4
        A
        B1
        B2
        B3
        B4

        %Constraint matrices
        %   E1*x + E2*u + E3*delta + E4 * z =< g5
        E1
        E2
        E3
        E4
        g5
        
        %Bounds
        %   bounds on continuous variables are necessary for mixed integer
        %   optimization, better bounds improve runtime and solutions
        xlb
        xub
        ulb
        uub
        zlb
        zub

        %State
        x      %current value of x
        xisbin %logical array signifying which states must be binary values

        %Variable names
        xname = {} %names of states
        uname = {} %names of states
        zname = {} %names of auxiliary variables
        dname = {} %names of delta variables

        %Number of variables
        nx %number of states
        ndelta %number of binary auxiliary variables
        nz %number of real auxiliary variables
        nu %number of inputs
        nconstraints %number of constraints
        
        %Other
        epsilon %machine precision value
        defopts %default optmizer options
        
    end
    methods
        %%
        %Constructor function
        %
        %no sanity checks are performed on input arguments
        %
        %required arguments: 
        %A,B1,B2,B3,B4,E1,E2,E3,E4,g5: matrices describing MLD model
        %xlb: lowerbound on state x
        %xub: upperbound on state x
        %zlb: lowerbound on slackvariables z
        %zub: upperbound on slackvariables z
        %nx: number of states
        %
        %optional arguments:
        %x0: starting state, default is zeros
        %xname: cellmatrix with names for states
        %uname: cellmatrix with names for inputs
        %dname: cellmatrix with names for binary auxiliary values
        %zname: cellmatrix with names for continuous auxiliary values
        %binvars: logical array with ones for binary states
        %epsilon: machine precision constant
        %opts: default optimizer options
        function obj = MLDModel(A,B1,B2,B3,B4,E1,E2,E3,E4,g5,xlb,xub,ulb,uub,zlb,zub,nx,varargin)
            
            if nargin < 18
                error("too few arguments")
            end
            
            obj.A = A;
            obj.B1 = B1;
            obj.B2 = B2;
            obj.B3 = B3;
            obj.B4 = B4;

            obj.E1 = E1;
            obj.E2 = E2;
            obj.E3 = E3;
            obj.E4 = E4;

            obj.g5 = g5;

            obj.xlb = xlb;
            obj.xub = xub;
            
            obj.ulb = ulb;
            obj.uub = uub;
            
            obj.zlb = zlb;
            obj.zub = zub;

            obj.nx = nx;
            obj.xisbin = logical([zeros(nx,1)]);
            
            obj.ndelta = size(B2,2);
            obj.nz = size(B3,2);
            obj.nu = size(B1,2);
            obj.nconstraints = size(g5,1);
            
            obj.x = zeros(obj.nx,1);
            obj.epsilon = 10^-6;
            obj.defopts = sdpsettings('verbose',0);
            
            obj.xname = {};
            for i = 1:obj.nx
                obj.xname = [obj.xname,"state"+string(i)]; 
            end
            
            obj.uname = {};
            for i = 1:obj.nu
                obj.uname = [obj.uname,"input"+string(i)]; 
            end
            
            obj.zname = {};
            for i = 1:obj.nz
                obj.zname = [obj.zname,"z"+string(i)]; 
            end
            
            obj.dname = {};
            for i = 1:obj.ndelta
                obj.dname = [obj.dname,"delta"+string(i)]; 
            end
    
            if (nargin > 18)
                i = 1;
                while(i <= size(varargin,2))
                    switch lower(varargin{i})
                        case "x0"; obj.x = varargin{i+1};
                        case "xname"; obj.xname = varargin{i+1};
                        case "uname"; obj.uname = varargin{i+1};
                        case "zname"; obj.zname = varargin{i+1};
                        case "dname"; obj.dname = varargin{i+1};
                        case "epsilon"; obj.epsilon = varargin{i+1};
                        case "opts"; obj.defopts = varargin{i+1};
                        case "binvars"; obj.xisbin = varargin{i+1};
                    end
                    i = i+2;
                end
            end
        end
        %%
        %Step method
        %
        %Advance mld model with a single step x = x_{k+1}
        function xkplus1 = step(obj,u,opts)
            if nargin<3
                opts = obj.defopts;
            end
            
            %Use yalmip to solve MILP
            xkp1 = sdpvar(obj.nx,1);

            if obj.ndelta ~= 0
                delta = binvar(obj.ndelta,1);
            else
                delta = 0;
            end
            if obj.nz ~= 0
                z = sdpvar(obj.nz,1);
            else
                z = 0;
            end
            
            C = [xkp1 == obj.A * obj.x + obj.B1 * u + obj.B2 * delta + obj.B3 * z + obj.B4]; %#ok<*NBRAK>
            C = [C, obj.E1 * obj.x + obj.E2* u + obj.E3 * delta + obj.E4 * z <= obj.g5];
            C = [C,obj.xlb <= xkp1 <= obj.xub];
            C = [C,obj.zlb <= z <= obj.zub]; 
            
            sol = optimize(C,[],opts);
            if sol.problem == 0
                obj.x = value(xkp1);
            else
                sol.info
                yalmiperror(sol.problem)
                error("Could not solve MILP problem")
            end
            xkplus1 = obj.x;
        end
        %%
        %Returns matrices that describe the N step ahead prediction
        %
        %matrices F1,F2,F3,M1,M2,M3 defined as
        %x(k+N) = M1*V + M2*x(k) + M3
        %F1*V <= F2 + F3*x(k)
        %
        %where
        %V = [U',DELTA',Z']'
        %U = [u(k)' u(k+1)' ... u(K+N-1)']'
        %DELTA = [delta(k)' delta(k+1)' ... delta(K+N-1)']'
        %Z = [z(k)' z(k+1)' ... z(K+N-1)']'
        function [F1,F2,F3,M1,M2,M3] = NstepMatrices(obj,N)
            [F1,F2,F3] = obj.NstepConstraintMatrices(N);
            [M1,M2,M3] = obj.NstepDynamicsMatrices(N);
        end
        %%
        %
        %
        %
        function [M1,M2,M3] = NstepDynamicsMatrices(obj,N)
            nx = obj.nx; %#ok<*PROPLC>
            nu = obj.nu;
            nd = obj.ndelta;
            nz = obj.nz;
            
            M11 = zeros(nx,nu*N);
            M11(:,end-nu+1:end) = obj.B1;
            M12 = zeros(nx,nd*N);
            M12(:,end-nd+1:end) = obj.B2;
            M13 = zeros(nx,nz*N);
            M13(:,end-nz+1:end) = obj.B3;
            
            M3 = zeros(nx,1);
            for i = 1:N-1
                M11(:,end-nu*(i+1)+1:end-nu*i) =  obj.A*M11(:,end-nu*(i)+1:end-nu*(i-1));
                M12(:,end-nd*(i+1)+1:end-nd*i) =  obj.A*M12(:,end-nd*(i)+1:end-nd*(i-1));
                M13(:,end-nz*(i+1)+1:end-nz*i) =  obj.A*M13(:,end-nz*(i)+1:end-nz*(i-1));
                
                M3 = M3 + obj.A^(i-1)*obj.B4;
            end
            
            M1 = [M11 M12 M13];
            M2 = obj.A^N;
            M3 = M3 + obj.A^(N-1)*obj.B4;
        end
        %%
        %
        %
        %
        function [F1,F2,F3] = NstepConstraintMatrices(obj,N)
            nx = obj.nx; %#ok<*PROPLC>
            ng = obj.nconstraints;

            F11 = kron(eye(N),obj.E2);
            F12 = kron(eye(N),obj.E3);
            F13 = kron(eye(N),obj.E4);

            F2 = kron(ones(N,1),obj.g5);
            
            F3 = zeros(N*ng,nx);
            F3(1:ng,:) = -obj.E1;
            
            vec = ones(N,1);
            for i = 1:N-1
                vi = i*ng;
                hi1 = 0;
                hi2 = 0;
                hi3 = 0;
                EA = obj.E1*obj.A^(i-1);
                F11i = EA*obj.B1;
                F12i = EA*obj.B2;
                F13i = EA*obj.B3;
                for j = 1:N-i
                    F11(1+vi:vi+ng,1+hi1:hi1+obj.nu) = F11i;
                    F12(1+vi:vi+ng,1+hi2:hi2+obj.ndelta) = F12i;
                    F13(1+vi:vi+ng,1+hi3:hi3+obj.nz) = F13i;
                    vi = vi+ng;
                    hi1 = hi1+obj.nu;
                    hi2 = hi2+obj.ndelta;
                    hi3 = hi3+obj.nz;
                end

                vec(i) = 0;
                F2 = F2 - kron(vec,obj.E1*obj.A^(i-1)*obj.B4);

                F3(i*ng+1:(i+1)*ng,:) = -obj.E1*obj.A^i;
            end
            
            F1 = [F11 F12 F13];
        end
        %%
        %Returns matrices that describe the state prediction over
        %the control horizon k..k+N
        %
        %%matrices M1,M2,M3 defined as
        %X = M1*V + M2*x(k) + M3
        %
        %where
        %X = [x(k+1)' x(k+2)' ... x(k+N)']'
        %V = [U',DELTA',Z']'
        %U = [u(k)' u(k+1)' ... u(K+N-1)']'
        %DELTA = [delta(k)' delta(k+1)' ... delta(K+N-1)']'
        %Z = [z(k)' z(k+1)' ... z(K+N-1)']'
        function [M1,M2,M3] = statePredictionMatrices(obj,N)
            M1 = zeros(obj.nx*N,(obj.nu+obj.ndelta+obj.nz)*N);
            M2 = zeros(obj.nx*N,obj.nx);
            M3 = zeros(obj.nx*N,1);
            
            
            
            nx = obj.nx;
            di = obj.nu*N;
            dz = (obj.nu+obj.ndelta)*N;
            for i = 1:N
                if i == 1
                    M1(1:nx,1:obj.nu) = obj.B1;
                    
                    M1(1:nx,di+1:di+obj.ndelta) = obj.B2;
                    
                    M1(1:nx,dz+1:dz+obj.nz) = obj.B3;
                    
                    M3(1:nx) = obj.B4;
                else
                    M1(1+(i-1)*nx:i*nx,1+obj.nu:di) = M1(1+(i-2)*nx:(i-1)*nx,1:obj.nu*(N-1));
                    M1(1+(i-1)*nx:i*nx,1:obj.nu) = obj.A*M1(1+(i-1)*nx:i*nx,1+obj.nu:2*obj.nu);

                    M1(1+(i-1)*nx:i*nx,1+di+obj.ndelta:di+obj.ndelta*N) = M1(1+(i-2)*nx:(i-1)*nx,1+di:di+obj.ndelta*(N-1));
                    M1(1+(i-1)*nx:i*nx,1+di:di+obj.ndelta) = obj.A*M1(1+(i-1)*nx:i*nx,1+di+obj.ndelta:di+2*obj.ndelta);

                    M1(1+(i-1)*nx:i*nx,1+dz+obj.nz:dz+obj.nz*N) = M1(1+(i-2)*nx:(i-1)*nx,1+dz:dz+obj.nz*(N-1));
                    M1(1+(i-1)*nx:i*nx,1+dz:dz+obj.nz) = obj.A*M1(1+(i-1)*nx:i*nx,1+dz+obj.nz:dz+2*obj.nz);

                    M3(1+(i-1)*nx:i*nx) = M3(1+(i-2)*nx:(i-1)*nx,:) + obj.A^i*obj.B4;
                end
                
                M2(1+(i-1)*nx:i*nx,:) = obj.A^i;
            end
        end
        %%
        %Add new constraint to the MLD model
        %
        %constraint must be of the form:
        %E1*x+E2*u+E3*delta+E4*z <= g5
        %
        %matrices can be empty
        function addConstraints(obj,E1,E2,E3,E4,g5,varargin)
            newconstraints = size(g5,1);
            obj.nconstraints = obj.nconstraints + newconstraints;
            
            if (nargin > 6)
                i = 1;
                while(i <= size(varargin,2))
                    switch lower(varargin{i})
                        case "xname"; xvars = varargin{i+1};
                        case "uname"; uvars = varargin{i+1};
                        case "zname"; zvars = varargin{i+1};
                        case "dname"; deltavars = varargin{i+1};    
                    end
                    i = i+2;
                end
            end
            
            if ~isempty(E1)
                if size(E1,2) == obj.nx
                    obj.E1 = [obj.E1;E1];
                else
                    if size(xvars,2) == size(E1,2)
                        tE1 = zeros(newconstraints,obj.nx);
                        n = 0;
                        for i = 1:obj.nx
                            for j = 1:size(xvars,2)
                                if strcmp(xvars{j},obj.xname{i})
                                    tE1(:,i) = E1(:,j);
                                    n = n+1;
                                    break
                                end
                            end
                        end
                        obj.E1 = [obj.E1;tE1];
                    end
                end
            else
                obj.E1 = [obj.E1;zeros(newconstraints,obj.nx)];
            end
            
            if ~isempty(E2)
                if size(E2,2) == obj.nu
                    obj.E2 = [obj.E2;E2];
                else
                    if size(uvars,2) == size(E2,2)
                        tE2 = zeros(newconstraints,obj.nu);
                        n = 0;
                        for i = 1:obj.nu
                            for j = 1:size(uvars,2)
                                if strcmp(uvars{j},obj.uname{i})
                                    tE2(:,i) = E2(:,j);
                                    n = n+1;
                                    break
                                end
                            end
                        end
                        obj.E2 = [obj.E2;tE2];
                    end
                end
            else
                obj.E2 = [obj.E2;zeros(newconstraints,obj.nu)];
            end
            
            if ~isempty(E3)
                if size(E3,2) == obj.ndelta
                    obj.E3 = [obj.E3;E3];
                else
                    if size(deltavars,2) == size(E3,2)
                        tE3 = zeros(newconstraints,obj.ndelta);
                        n = 0;
                        for i = 1:obj.ndelta
                            for j = 1:size(deltavars,2)
                                if strcmp(deltavars{j},obj.dname{i})
                                    tE3(:,i) = E3(:,j);
                                    n = n+1;
                                    break
                                end
                            end
                        end
                        obj.E3 = [obj.E3;tE3];
                    end
                end
            else
                obj.E3 = [obj.E3;zeros(newconstraints,obj.ndelta)];
            end
            
            if ~isempty(E4)
                if size(E4,2) == obj.nz
                    obj.E4 = [obj.E4;E4];
                else
                    if size(zvars,2) == size(E4,2)
                        tE4 = zeros(newconstraints,obj.nz);
                        n = 0;
                        for i = 1:obj.nz
                            for j = 1:size(zvars,2)
                                if strcmp(zvars{j},obj.zname{i})
                                    tE4(:,i) = E4(:,j);
                                    n = n+1;
                                    break
                                end
                            end
                        end
                        obj.E4 = [obj.E4;tE4];
                    end
                end
            else
                obj.E4 = [obj.E4;zeros(newconstraints,obj.nz)];
            end
            
           obj.g5 = [obj.g5;g5];
        end
        %%
        %
        %
        %
        function addStates(obj,nx,lb,ub,A1,A2,B1,B2,B3,B4,E1,varnames)
            if isempty(A1)
                A1 = zeros(obj.nx,nx);
            end
            if isempty(A2)
                A2 = zeros(nx,obj.nx+nx);
            end
            obj.A = [obj.A A1;A2];
            
            if isempty(B1)
                B1 = zeros(nx,obj.nu);
            end
            obj.B1 = [obj.B1;B1];
            
            if isempty(B2)
                B2 = zeros(nx,obj.ndelta);
            end
            obj.B2 = [obj.B2;B2];

            if isempty(B3)
                B3 = zeros(nx,obj.nz);
            end
            obj.B3 = [obj.B3;B3];

            if isempty(B4)
               B4 = zeros(nx,1);
            end
            obj.B4 = [obj.B4;B4];
            
            if isempty(E1)
                E1 = zeros(obj.nconstraints,nx);
            end
            obj.E1 = [obj.E1,E1];
            
            if (nargin > 11)
                if size(varnames) ~= nx
                    error("size of names list is wrong")
                end
                obj.xname = [obj.xname,varnames];
            else
                for i = 1:nx
                    obj.xname = [obj.xname,{"state"+string(obj.nx+i)}]; 
                end
            end
            obj.nx = obj.nx+nx;
            obj.xlb = [obj.xlb;lb];
            obj.xub = [obj.xub;ub];
        end
        %%
        %
        %
        %
        function addInputs(obj,nu,lb,ub,varnames)
            obj.B1 = [obj.B1,zeros(obj.nx,nu)];
            obj.E2 = [obj.E2,zeros(obj.nconstraints,nu)];
            
            if (nargin > 4)
                if size(varnames) ~= nu
                    error("size of names list is wrong")
                end
                obj.uname = [obj.uname,varnames];
            else
                for i = 1:nu
                    obj.uname = [obj.uname,{"input"+string(obj.nu+i)}]; 
                end
            end
            obj.nu = obj.nu+nu;
            obj.ulb = [obj.ulb;lb];
            obj.uub = [obj.uub;ub];
        end
        %%
        %
        %
        %
        function addDeltas(obj,nd,varnames)
            obj.B2 = [obj.B2,zeros(obj.nx,nd)];
            obj.E3 = [obj.E3,zeros(obj.nconstraints,nd)];
            
            if (nargin > 2)
                if size(varnames,2) ~= nd
                    error("size of names list is wrong")
                end
                obj.dname = [obj.dname,varnames];
            else
                for i = 1:nd
                    obj.dname = [obj.dname,{"delta"+string(obj.nd+i)}]; 
                end
            end
            obj.ndelta = obj.ndelta+nd;
        end
        %%
        %
        %
        %
        function addAuxvariables(obj,nz,lb,ub,varnames)
            obj.B3 = [obj.B3,zeros(obj.nx,nz)];
            obj.E4 = [obj.E4,zeros(obj.nconstraints,nz)];
            obj.zlb = [obj.zlb;lb];
            obj.zub = [obj.zub;ub];
            
            if (nargin > 4)
                if size(varnames,2) ~= nz
                    error("size of names list is wrong")
                end
                obj.zname = [obj.zname,varnames];
            else
                for i = 1:nz
                    obj.zname = [obj.zname,{"z"+string(obj.nz+i)}]; 
                end
            end
            obj.nz = obj.nz+nz;
        end
        %%
        %Adds binary equivalence constraint
        %F1*x+F2*u+F3*z+f4 <= 0 <-> [delta = 1]
        %
        function addBinaryEquivalence(obj,F1,F2,F3,f4,dname,varargin)
            xvars = [];
            uvars = [];
            zvars = [];
            if (nargin > 6)
                i = 1;
                while(i <= size(varargin,2))
                    switch lower(varargin{i})
                        case "xname"; xvars = varargin{i+1};
                        case "uname"; uvars = varargin{i+1};
                        case "zname"; zvars = varargin{i+1};
                    end
                    i = i+2;
                end
            end
            
            if F1 == 0
                F1 = [];
            end
            if F2 == 0
                F2 = [];
            end
            if F3 == 0
                F3 = [];
            end
            
            if (~isempty(F1)) && (isempty(xvars))
                xvars = obj.xname;
            end
            if (~isempty(F2)) && (isempty(uvars))
                uvars = obj.uname;
            end
            if (~isempty(F3)) && (isempty(zvars))
                zvars = obj.zname;
            end
            
            lb = obj.getLowerBounds([xvars,uvars,zvars]);
            ub = obj.getUpperBounds([xvars,uvars,zvars]);
            [~,m] = linprog([F1,F2,F3]',[],[],[],[],lb,ub);
            [~,M] = linprog([-F1,-F2,-F3]',[],[],[],[],lb,ub);
            m = m+f4;
            M = -M+f4;
            
            obj.addConstraints(F1,F2,M,F3,M,"xname",xvars,"uname",uvars,"zname",zvars,"dname",dname);
            obj.addConstraints(-F1,-F2,(m-obj.epsilon),-F3,-obj.epsilon,"xname",xvars,"uname",uvars,"zname",zvars,"dname",dname);
        end
        %%
        %adds new auxiliary variable z with constraints so
        %z = delta*(F1*x+F2*u+F3*z+f4)
        %
        function addBinaryContinuousMultiplication(obj,F1,F2,F3,f4,dname,varargin)
            xvars = [];
            uvars = [];
            zvars = [];
            name = [];
            if (nargin > 6)
                i = 1;
                while(i <= size(varargin,2))
                    switch lower(varargin{i})
                        case "xname"; xvars = varargin{i+1};
                        case "uname"; uvars = varargin{i+1};
                        case "zname"; zvars = varargin{i+1};
                        case "variablename"; name = varargin{i+1};
                    end
                    i = i+2;
                end
            end
            
            if F1 == 0
                F1 = [];
            end
            if F2 == 0
                F2 = [];
            end
            if F3 == 0
                F3 = [];
            end
            if isempty(f4)
                f4 = 0;
            end
            
            if (~isempty(F1)) && (isempty(xvars))
                xvars = obj.xname;
            end
            if (~isempty(F2)) && (isempty(uvars))
                uvars = obj.uname;
            end
            if (~isempty(F3)) && (isempty(zvars))
                zvars = obj.zname;
            end
            
            lb = obj.getLowerBounds([xvars,uvars,zvars]);
            ub = obj.getUpperBounds([xvars,uvars,zvars]);
            [~,m] = linprog([F1,F2,F3]',[],[],[],[],lb,ub);
            [~,M] = linprog([-F1,-F2,-F3]',[],[],[],[],lb,ub);
            m = m+f4;
            M = -M+f4;
            
            if m > 0
                m = 0;
            end
            if M < 0
                M = 0;
            end
            
            if ~isempty(name)
                obj.addAuxvariables(1,m,M,name);
            else
                obj.addAuxvariables(1,m,M);
            end
            
            obj.addConstraints([],[],[-M;m],[1;-1],[0;0],"dname",dname,"zname",string(obj.zname{end}));
            obj.addConstraints(-F1,-F2,-m,[-F3,1],f4-m,"xname",xvars,"uname",uvars,"dname",dname,"zname",[zvars,string(obj.zname{end})]);
            obj.addConstraints(F1,F2,M,[F3,-1],-f4+M,"xname",xvars,"uname",uvars,"dname",dname,"zname",[zvars,string(obj.zname{end})]);
        end
        
        %%
        %
        %
        function lb = getLowerBounds(obj,varnames)
            lb = [];
            for i = 1:length(varnames)
                xi = find(obj.xname == varnames{i},1);
                if isempty(xi)
                    ui = find(obj.uname == varnames{i},1);
                    if isempty(ui)
                        zi = find(obj.zname == varnames{i},1);
                        if isempty(zi)
                            error("incorrect variable name");
                        else
                            lb = [lb;obj.zlb(zi)];
                        end
                    else
                        lb = [lb;obj.ulb(ui)];
                    end
                else
                    lb = [lb;obj.xlb(xi)];
                end
            end
        end
        %%
        %
        %
        function ub = getUpperBounds(obj,varnames)
            ub = [];
            for i = 1:length(varnames)
                xi = find(obj.xname == varnames{i},1);
                if isempty(xi)
                    ui = find(obj.uname == varnames{i},1);
                    if isempty(ui)
                        zi = find(obj.zname == varnames{i},1);
                        if isempty(zi)
                            error("incorrect variable name");
                        else
                            ub = [ub;obj.zub(zi)];
                        end
                    else
                        ub = [ub;obj.uub(ui)];
                    end
                else
                    ub = [ub;obj.xub(xi)];
                end
            end
        end
    end
    methods(Static)
        %%
        %
        %
        function model = new(varargin)
            model = MLDModel([],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],0,varargin);
        end
    end
end