classdef microgridTreePredictor
    
    properties
        treeClassifiers
        N
        ndelta
        feasibilityfunction
    end
    
    methods
        function obj = microgridTreePredictor(classifiers,N,ndelta,feasibilityfunction)
            obj.treeClassifiers = classifiers;
            obj.N = N;
            obj.ndelta = ndelta;
            obj.feasibilityfunction = feasibilityfunction;
        end
        
        function deltas = predict(obj,z)
            deltas = zeros(obj.N,obj.ndelta);
            zk = cellfun(@transpose,z,'UniformOutput',false); 
            for i = 1:(obj.N-1)
                deltas(i,:) = str2num(cell2mat(obj.treeClassifiers{i}.predict(cell2mat(zk))));
                [~,deltas(i,:)] = obj.feasibilityfunction(deltas(i,:),z{5}(i),z{6}(i));
            end
            deltas = reshape(deltas',[obj.ndelta*obj.N,1]);
        end
    end
end

