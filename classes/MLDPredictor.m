classdef (Abstract) MLDPredictor < handle
    properties
    end
    methods (Abstract)
        deltas = predict(obj,z,feasibility)
    end
end