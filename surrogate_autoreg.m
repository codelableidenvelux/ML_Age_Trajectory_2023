%% C. Autoregressive model
% This is built following the equations in
% https://www-sciencedirect-com.ezproxy.leidenuniv.nl/science/article/pii/016727899290102S?source=UBLbookmarklet

function ys = surrogate_autoreg(x, n)
    mu = mean(x);
    v = var(x);
    ys = zeros(n, length(x));
    ac = autocorr(x);
    a1 = ac(2);
    a0 = mu * (1 - a1);
    sigma = sqrt(v * (1 - a1^2));
    for i = 1:n

        ys(i, :) = my_ar(x(1), a0, a1, sigma, length(x));
    end
end

function y = my_ar(x0, a0, a1, sigma, n)

    xt = @(xt_1) a0 + a1 * xt_1 + randn() * sigma;
    
    y = zeros(1,n);
    y(1) = x0;
    
    for i = 2:n
    y(i) = xt(y(i - 1));
    end
end