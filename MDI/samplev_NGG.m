function [v] = samplev_NGG(lambda, sumv, n, K, M, a)

%maxstar = (lambda + sumv)^(K * a - (n - 1)) * exp(- M / a * (lambda + sumv)^a);
logmaxstar = (K * a - (n - 1)) * log(lambda + sumv) - M / a * (lambda + sumv)^a;

u1 = rand;
logu = log(u1) + logmaxstar;

if ( n > 5 )
    z = M + n - 1 + K * a;
    start = log(u1^(-1/z) - 1) + log(lambda + sumv);
else
    start = 0;
end

if ( (isnan(start) == 0) || (isinf(start) == 0 ) )
    start = 0;
end

%[logv, fval, exitflag] = fzero('f1_NGG', start, optimset('Display','iter'), lambda, sumv, n, K, M, a, logu);
[logv, fval, exitflag] = fzero('f1_NGG', start, [], lambda, sumv, n, K, M, a, logu);

if ( exitflag ~= 1 )
    lambda
    sumv
    n
    K
    M
    a
    u
    logv
    fval
end
    

v = exp(logv);

