classdef microgridHeuristicPredictor
    
    properties
        N
        mgrid
    end
    
    methods
        function obj = microgridHeuristicPredictor(N,mgrid)
            obj.N = N;
            obj.mgrid = mgrid;
        end
        
        function [deltas, dcase] = predict(obj,z)
            deltas = zeros(obj.N*obj.mgrid.ndelta,1);
            csale = z{3};
            cpur = z{2};
            cprod = z{4};
            Pload = z{5};
            Pres = z{6};
            dcase = zeros(obj.N,1);
            for i = 1:obj.N
                r = 1+(i-1)*4:4+(i-1)*4;
                dcase(i) = 0;
                c1 = Pload(i) <= Pres(i);
                c2 = (obj.mgrid.uub(2)+obj.mgrid.uub(3) >= Pload(i)-Pres(i)) & (Pload(i) > Pres(i));
                c3 = obj.mgrid.uub(2)+obj.mgrid.uub(3) < Pload(i)-Pres(i);
                c4 = (cprod(i) < csale(i)) & (csale(i) <= cpur(i));
                c5 = (csale(i) < cprod(i)) & (cprod(i) <= cpur(i));
                c6 = (csale(i) < cpur(i)) & (cpur(i) <= cprod(i));
                if c1
                    if c4
                        dcase(i) = 1;
                    else
                        dcase(i) = 4;
                    end
                elseif c2
                    if c6
                        dcase(i) = 2;
                    else
                        dcase(i) = 5;
                    end
                else
                    if c6
                        dcase(i) = 5;
                    else
                        dcase(i) = 3;
                    end
                end
                switch dcase(i)
                    case 1
                        deltas(r) = [0 1 1 1];
                    case 2
                        deltas(r) = [0 1 1 1];
                    case 3
                        deltas(r) = [1 1 1 0];
                    case 4
                        deltas(r) = [0 0 0 1];
                    case 5
                        deltas(r) = [0 1 1 1];
                end
            end
        end
    end
end

