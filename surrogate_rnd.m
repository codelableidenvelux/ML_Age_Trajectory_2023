%% A. randomization

function ys = surrogate_rnd(x, n)

    x1 = x(randperm(length(x)));
    ys = zeros(n, length(x));
    
    for i = 1:n
        ys(i, :) = x1(randi(length(x1), 1, length(x1)));
    end
end