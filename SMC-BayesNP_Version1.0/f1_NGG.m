function [logf1] = f1_NGG(logv, lambda, sumv, n, K, M, a, logu)

v = exp(logv);

%f1 = (lambda + sumv + v)^(K * a - (n - 1)) * exp(- M / a * (lambda + sumv + v)^a) - u;

logf1 = (K * a - (n - 1)) * log(lambda + sumv + v) - M / a * (lambda + sumv + v)^a - logu;
