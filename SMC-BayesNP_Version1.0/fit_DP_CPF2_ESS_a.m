function [holdnumbclust, holda] = fit_DP_CPF2_ESS_a(data, mu, sigmasq, M, burnin, numbofits, every, numbofparts, adap)

a = 0.1;
aprior = [1 1];

s = 1:length(data);
sumv = zeros(1, length(data));

lambda = gamrnd(M, 1);
sumv(1) = sumv(1) - log(rand) / lambda;

for i = 2:length(data)
    lambda = gamrnd(M + i - 1, 1 ./ (1 + sumv(i - 1)));
    sumv(i) = sumv(i - 1) - log(rand) / lambda;
end

holdnumbclust = zeros(1, numbofits);
holda = zeros(1, numbofits);

aaccept = 0;
acount = 0;
logasd = log(0.01);

for it = 1:(burnin + numbofits * every)
    
    [s, sumv] = filter1_DP_conj_CPF2_ESS(s, sumv, data, mu, sigmasq, a, M, numbofparts, adap);
    
    trans = log(a) - log(1 - a);
    newtrans = trans + exp(logasd) * randn;
    newa = exp(newtrans) / (1 + exp(newtrans));
     
    loglike = - 0.5 * length(data) * log(a) - 0.5 * length(unique(s)) * log(1 - a);
    loglike = loglike - 0.5 * sum(data.^2) / (a * sigmasq) - 0.5 * length(unique(s)) * mu^2 / ((1 - a) * sigmasq);
    for i = 1:length(unique(s))
        loglike = loglike - 0.5 * log(sum(s==i) / a + 1 / (1 - a));
        loglike = loglike + 0.5 / sigmasq * (sum(data(s==i)) / a + mu / (1 - a))^2 / (sum(s==i) / a + 1 / (1 - a));       
    end
    
    newloglike = - 0.5 * length(data) * log(newa) - 0.5 * length(unique(s)) * log(1 - newa);
    newloglike = newloglike - 0.5 * sum(data.^2) / (newa * sigmasq) - 0.5 * length(unique(s)) * mu^2 / ((1 - newa) * sigmasq);
    for i = 1:length(unique(s))
        newloglike = newloglike - 0.5 * log(sum(s==i) / newa + 1 / (1 - newa));
        newloglike = newloglike + 0.5 / sigmasq * (sum(data(s==i)) / newa + mu / (1 - newa))^2 / (sum(s==i) / newa + 1 / (1 - newa));       
    end
    
    logaccept = newloglike - loglike - log(1 / newa + 1 / (1 - newa)) + log(1 / a + 1 / (1 - a));
    logaccept = logaccept + (aprior(1) - 1) * (log(newa) - log(a)) + (aprior(2) - 1) * (log(1 - newa) - log(1 - a));
    
    accept = 1;
    if ( logaccept < 0 )
        accept = exp(logaccept);
    end
    
    if ( rand < accept )
        a = newa;
    end

    logasd = logasd + 1 / it^(0.55) * (accept - 0.234);
    aaccept = aaccept + accept;
    acount = acount + 1;
    
    if ( (it > burnin) && (mod(it - burnin, every) == 0) )
        holdnumbclust((it - burnin) / every) = length(unique(s));
        holda((it - burnin) / every) = a;
    end
end

