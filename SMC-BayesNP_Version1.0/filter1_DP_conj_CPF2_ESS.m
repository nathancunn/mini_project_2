function [s, sumv] = filter1_DP_conj_CPF2_ESS(s_old, sumv_old, data, mu, sigmasq, a, M, numbofparts, adap)

n = length(data);

sumy = cell(1, numbofparts);
sumysq = cell(1, numbofparts);
nj = cell(1, numbofparts);
sumv = zeros(numbofparts, n);
s = zeros(numbofparts, n);


logweight = zeros(1, numbofparts);
indx = zeros(1, numbofparts);

for it = 1:numbofparts
    sumy{1, it} = data(1);
    sumysq{1, it} = data(1)^2;
    nj{1, it} = 1;
    logweight(it) = 0;
    s(it, 1) = 1;
end

lambda = gamrnd(M, 1, numbofparts - 1, 1);
sumv(1:(end-1), 1) = sumv(1:(end-1), 1) - log(rand(numbofparts - 1, 1)) ./ lambda;
sumv(numbofparts, 1) = sumv_old(1);

for i = 2:n
    weight = exp(logweight - max(logweight));
    weight = weight / sum(weight);

    if ( adap == 1 )
        ESS = 1 / sum(weight.^2);
        if ( ESS < 0.5 * numbofparts )
            update = 1;
        else
            update = 0;
        end
    else
        update = 1;
    end
        
    if ( update == 1 )
        fprob = cumsum(weight);
        fprob = fprob / fprob(end);
        
        for part = 1:(numbofparts-1)
            u = rand;
            indx(part) = 1;
            while ( u > fprob(indx(part)) )
                indx(part) = indx(part) + 1;
            end
        end
        indx(numbofparts) = numbofparts;
        
        s(:, 1:(i-1)) = s(indx, 1:(i-1));
        nj = nj(1, indx);
        sumy = sumy(1, indx);
        sumysq = sumysq(1, indx);
        sumv = sumv(indx, 1:(i-1));
        
        logweight = zeros(1, numbofparts);
    end
    
    lambda = gamrnd(M + i - 1, 1 ./ (1 + sumv(1:(end-1), i - 1)));
    sumv(1:(end-1), i) = sumv(1:(end-1), i - 1) - log(rand(numbofparts - 1, 1)) ./ lambda;
    for it = 1:(numbofparts-1)
        prob = [nj{1, it} M] / (sum(nj{1, it}) + M);
        mustar = [(sumy{1, it} / a + mu / (1 - a)) ./ (nj{1, it} / a + 1 / (1 - a)) mu];
        varstar = [sigmasq ./ (nj{1, it} / a + 1 / (1 - a)) sigmasq*(1-a)];
        logprob = - 0.5 * (data(i) - mustar).^2 ./ (a * sigmasq + varstar) - 0.5 * log(a * sigmasq + varstar);
        fprob = cumsum(prob .* exp(logprob - max(logprob)));
        logweight(it) = logweight(it) + log(fprob(end)) + max(logprob);
        fprob = fprob / fprob(end);
        u1 = rand;
        s(it, i) = 1;
        while ( fprob(s(it, i)) < u1 )
            s(it, i) = s(it, i) + 1;
        end
    end
    sumv(numbofparts, i) = sumv_old(i);
    s(numbofparts, i) = s_old(i);
    prob = [nj{1, numbofparts} M] / (sum(nj{1, numbofparts}) + M);
    mustar = [(sumy{1, numbofparts} / a + mu / (1 - a)) ./ (nj{1, numbofparts} / a + 1 / (1 - a)) mu];
    varstar = [sigmasq ./ (nj{1, numbofparts} / a + 1 / (1 - a)) sigmasq*(1-a)];
    logprob = - 0.5 * (data(i) - mustar).^2 ./ (a * sigmasq + varstar) - 0.5 * log(a * sigmasq + varstar);
    fprob = cumsum(prob .* exp(logprob - max(logprob)));
    logweight(numbofparts) = logweight(numbofparts) + log(fprob(end)) + max(logprob);
    
        
    for it = 1:numbofparts
        if ( s(it, i) > length(nj{1, it}) )
            nj{1, it}(s(it, i)) = 1;
            sumy{1, it}(s(it, i)) = data(i);
            sumysq{1, it}(s(it, i)) = data(i)^2;
        else
            nj{1, it}(s(it, i)) = nj{1, it}(s(it, i)) + 1;
            sumy{1, it}(s(it, i)) = sumy{1, it}(s(it, i)) + data(i);
            sumysq{1, it}(s(it, i)) = sumysq{1, it}(s(it, i)) + data(i)^2;
        end
    end
end

fprob = cumsum(exp(logweight  - max(logweight)));
fprob = fprob / fprob(end);

u = rand;
counter = 1;
while ( u > fprob(counter) )
    counter = counter + 1;
end
s = s(counter, :);
sumv = sumv(counter, :);

