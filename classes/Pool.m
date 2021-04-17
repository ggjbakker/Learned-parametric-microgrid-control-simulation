classdef Pool < handle
    properties
        genes
        data
        
        constants
        nconstants
        nvariables
        nsamples
        nhorizon
        nchromosomes
        
        genlength
        taillength
        totalgenlength
        populationsize
       
        f
    end
    
    methods
        function obj = Pool(populationsize,gensize,nchroms,data,constants)
            obj.populationsize = populationsize;
            obj.taillength = gensize+1;
            obj.genlength = gensize+obj.taillength;
            obj.totalgenlength = obj.genlength*nchroms*3;
            obj.nchromosomes = nchroms;
            obj.data = data;
            obj.constants = constants;
            obj.nconstants = length(constants);
            obj.nsamples = size(data.v,3);
            obj.nhorizon = size(data.v,2);
            obj.nvariables = size(data.v,1);
            obj.f = zeros(obj.nhorizon,obj.nchromosomes);
        end
        function fit = testFitness(obj,gen)
            fit = 0;
            expressions = cell(3,1);
            for m = 0:obj.nchromosomes-1
                chrom = gen((m*obj.genlength)*3+1:3*(m+1)*obj.genlength);
                expressions{m+1} = Tree.fromGen(char(chrom));
            end
            for i = 1:obj.nsamples
                v = obj.data.v(:,:,i);
                for m = 1:obj.nchromosomes
                    for k = 1:obj.nhorizon
                        obj.f(k,m) = expressions{m}.evaluate(obj.constants,v(:,k));
                    end
                end
                %uhat = obj.f*(obj.f\obj.data.u(:,i));
                %fit = fit + sumsqr(obj.data.u(:,i)-uhat)/obj.nsamples;
                annih = eye(obj.nhorizon)-obj.f*pinv(obj.f);
                se = obj.data.u(:,i)'*annih*obj.data.u(:,i);
                fit = fit + se/obj.nsamples;
            end
            fit = fit/(obj.nhorizon-obj.nchromosomes);
        end
        function populate(obj)
            for i = 1:obj.populationsize
                obj.genes(:,i) = obj.randomGen();
            end
        end
        function [meanfitness,best] = step(obj)
            fitness = obj.evaluatePopulation();
            meanfitness = mean(fitness(:,2));
            best = min(fitness(:,2));
            successors = obj.selectSuccesors(fitness);
            obj.genes = obj.adaptPopulation(successors);
        end
        function fitness = evaluatePopulation(obj)
            fitness = zeros(obj.populationsize,1);
            parfor i = 1:obj.populationsize
                fitness(i) = obj.testFitness(char(obj.genes(:,i)));
            end
            fitness = [[1:obj.populationsize]',fitness];
        end
        function successors = selectSuccesors(obj,fitness)
            successors = zeros(obj.populationsize/4,1);
            [~,best] = min(fitness(:,2));
            successors(1) = best;
            fitness(best,:) = [];
            fitness(:,2) = fitness(:,2)/mean(fitness(:,2));
            fitness(:,2) = 1./fitness(:,2);
            fitness = fitness(~any( isinf( fitness ), 2 ),: );
            for i = 2:(obj.populationsize/2)
                probabilities = fitness(:,2)./sum(fitness(:,2));
                r = rand;
                for k = 1:length(probabilities)
                    if r < probabilities(k)
                        successors(i) = fitness(k,1);
                        fitness(k,:) = [];
                        break
                    else
                        r = r-probabilities(k);
                    end
                end
            end
        end
        function newgenes = adaptPopulation(obj,successors)
            pmutate = 0.7;
            ns = length(successors);
            newgenes = repmat(' ',obj.totalgenlength,ns);
            for i = 1:ns
                newgenes(:,i) = obj.genes(:,successors(i));
            end
            for i = 1:(obj.populationsize-ns)
                if rand < pmutate
                    r = randi([1 ns]);
                    newgenes(:,ns+i) = obj.mutate(obj.genes(:,r));
                else
                    r = randi([1 ns],1,2);
                    while r(1) == r(2)
                        r = randi([1 ns],1,2);
                    end
                    newgenes(:,ns+i) = obj.crossover(obj.genes(:,r(1)),obj.genes(:,r(2)));
                end
            end
        end
        function gen = mutate(obj,gen)
            r = randi([1 obj.genlength*obj.nchromosomes]);
            if 0 < mod(r,obj.genlength) && mod(r,obj.genlength) < (obj.genlength-obj.taillength)
                gen((r-1)*3+1:r*3) = obj.randomLetter(1);
            else
                gen((r-1)*3+1:r*3) = obj.randomLetter(0);
            end
        end
        function gen = crossover(obj,gen1,gen2)
            r = randi([1 (obj.totalgenlength-3)/3]);
            gen = [gen1(1:r*3);gen2(r*3+1:end)];
        end
        function gen = randomGen(obj)
            gen = blanks(obj.totalgenlength);
            for m = 0:obj.nchromosomes-1
                i = m*obj.genlength*3;
                for p = 0:obj.genlength-1
                    if p < (obj.genlength-obj.taillength)
                        gen(p*3+i+1:3*(p+1)+i) = obj.randomLetter(1);
                    else
                        gen(p*3+i+1:3*(p+1)+i) = obj.randomLetter(0);
                    end
                end
            end
        end
        function letter = randomLetter(obj,inHead)
            Pfn = 0.6;
            PI = 0;
            if inHead && (rand < Pfn)
                p = rand;
                if p < 0.3
                    letter = 'f+_';
                elseif p < 0.6
                    letter = 'f-_';
                elseif p < 0.8
                    letter = 'f*_';
                else
                    letter = 'f+_';
                end
            else
                p = rand;
                if p < (1-PI)
                    tot = obj.nconstants+obj.nvariables;
                    number = floor(rand*tot)+1;
                    if number <= obj.nconstants
                        letter = ['c',char(num2str(number,'%.2i'))];
                    else
                        number = number - obj.nconstants;
                        letter = ['v',char(num2str(number,'%.2i'))];
                    end
                else
                    number = floor(rand*5)+1;
                    letter = ['I',char(num2str(number,'%.2i'))];
                end
            end
        end
    end
end

