function [s] = algorithm1_DP(data, mu, sigmasq, a, M, numbofparts)

n = length(data);

sumy = cell(1, numbofparts);
sumysq = cell(1, numbofparts);
nj = cell(1, numbofparts);
sumv = zeros(numbofparts, 1);
s = zeros(numbofparts, n);

logweight = zeros(1, numbofparts);
partstar = zeros(1, numbofparts);

for i = 1:n
    disp(['i = ' num2str(i)]);
    
    for it = 1:numbofparts
        
        if ( i == 1 )
            sumy{1, it} = data(i);
            sumysq{1, it} = data(i)^2;
            nj{1, it} = 1;
            logweight(it) = 0;
            s(it, i) = 1;
        else
            prob = [nj{1, it} M] / (sum(nj{1, it}) + M);
            mustar = [(mu / (1 - a) + sumy{1, it} / a) ./ (nj{1, it} / a + 1 / (1 - a)) mu];
            varstar = [sigmasq ./ (nj{1, it} / a + 1 / (1 - a)) sigmasq*(1-a)];
            logprob = - 0.5 * (data(i) - mustar).^2 ./ (a * sigmasq + varstar) - 0.5 * log(a * sigmasq + varstar);
            fprob = cumsum(prob .* exp(logprob - max(logprob)));
            logweight(it) = log(fprob(end)) + max(logprob);
            fprob = fprob / fprob(end);
            u1 = rand;
            s(it, i) = 1;
            while ( fprob(s(it, i)) < u1 )
                s(it, i) = s(it, i) + 1;
            end
            if ( s(it, i) == length(nj{1, it}) + 1 )
                nj{1, it} = [nj{1, it} 1];
                sumy{1, it} = [sumy{1, it} data(i)];
                sumysq{1, it} = [sumysq{1, it} data(i)^2];
            else
                nj{1, it}(s(it, i)) = nj{1, it}(s(it, i)) + 1;
                sumy{1, it}(s(it, i)) = sumy{1, it}(s(it, i)) + data(i);
                sumysq{1, it}(s(it, i)) = sumysq{1, it}(s(it, i)) + data(i)^2;
            end
        end
    end
    
    fprob = cumsum(exp(logweight - max(logweight)));
    fprob = fprob / fprob(end);
        
    u1 = rand / numbofparts;
    j = 1;
    it = 1;
    while ( j <= numbofparts )
        while ( (u1 < fprob(j)) && (j <= numbofparts) )
            partstar(it) = j;
            u1 = u1 + 1 / numbofparts;
            it = it + 1;
        end
        j = j + 1;
    end
    
    sumy = sumy(1, partstar);
    sumysq = sumysq(1, partstar);
    nj = nj(1, partstar);
    sumv = sumv(partstar);
    s = s(partstar, :);
end
