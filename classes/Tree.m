classdef Tree < handle
    properties
        left
        right
        rep
        isTerminal
        hasChildren
        value
    end
    
    methods
        function obj = Tree(rep)
            obj.rep = rep;
            if obj.rep(1) == 'f'
                obj.isTerminal = 0;
                obj.hasChildren = 0;
            else
                obj.isTerminal = 1;
                obj.value = str2double(obj.rep(2:3));
            end
        end
        
        function  result = evaluate(obj,c,v)
            if obj.rep(1) == 'f'
                if obj.rep(2) == '+'
                    result = obj.left.evaluate(c,v)+obj.right.evaluate(c,v);
                elseif obj.rep(2) == '-'
                    result = obj.left.evaluate(c,v)-obj.right.evaluate(c,v);
                elseif obj.rep(2) == '*'
                    result = obj.left.evaluate(c,v)*obj.right.evaluate(c,v);
                else
                    result = obj.left.evaluate(c,v)/obj.right.evaluate(c,v);
                end
            elseif obj.rep(1) == 'c'
                result = c(obj.value);
            elseif obj.rep(1) == 'v'
                result = v(obj.value);
            else
                result = obj.value;
            end
        end
        
        function string = representation(obj)
            if obj.rep(1) == 'f'
                string = obj.rep + obj.left.representation() + obj.right.representation();
            else
                string = obj.rep;
            end 
        end
        
        function gen = expand(obj,gen) 
            if not(obj.isTerminal)
                if obj.hasChildren
                    gen = obj.left.expand(gen);
                    gen = obj.right.expand(gen);
                    if obj.left.isTerminal && obj.right.isTerminal
                        obj.isTerminal = 1;
                    end
                else
                    obj.left = Tree(gen(1:3));
                    obj.right = Tree(gen(4:6));
                    obj.hasChildren = 1;
                    if obj.left.isTerminal && obj.right.isTerminal
                        obj.isTerminal = 1;
                    end
                    gen = gen(7:end);
                end
            end
        end
    end
    methods(Static)
        function [tree, gen] = fromGen(gen)
            tree = Tree(gen(1:3));
            gen = gen(4:end);
            while not(tree.isTerminal)
                gen = tree.expand(gen);
            end
        end
   end
end

