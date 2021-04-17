function [values] = evaluateGen(gen,r,c,v)
    reps = gen(1:r*3);
    gen = gen(r*3+1:end);
    rnext = sum(reps(1:3:end)=='f')*2;
    if rnext ~= 0
        nextlayer = evaluateGen(gen,rnext,c,v);
    end
    values = zeros(r,1);
    k = 1;
    for i = 1:r
        rep = reps(i*3-2:i*3);
        if rep(1) == 'f'
            if rep(2) == '+'
                values(i) = nextlayer(k) + nextlayer(k+1);
            elseif rep(2) == '-'
                values(i) = nextlayer(k) - nextlayer(k+1);
            elseif rep(2) == '*'
                values(i) = nextlayer(k) * nextlayer(k+1);
            else
                values(i) = nextlayer(k) / nextlayer(k+1);
            end
            k = k+2;
        elseif rep(1) == 'c'
            values(i) = c(round(str2double(rep(2:3))));
        elseif rep(1) == 'v'
            values(i) = v(round(str2double(rep(2:3))));
        else
            values(i) = round(str2double(rep(2:3)));
        end
    end
end

