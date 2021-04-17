classdef (Abstract) MLDController < handle
    properties
    end
    methods (Abstract)
        u = control(obj,x0)
    end
end