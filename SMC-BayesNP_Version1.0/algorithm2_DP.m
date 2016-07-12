function [s] = algorithm2_DP(data, mu, sigmasq, a, M, m, numbofparts)

n = length(data);

sumy = cell(1, numbofparts);
sumysq = cell(1, numbofparts);
nj = cell(1, numbofparts);
s = zeros(numbofparts, n);
theta = cell(1, numbofparts);

logweight = zeros(1, numbofparts);
partstar = zeros(1, numbofparts);

for i = 1:n
    disp(['i = ' num2str(i)]);
    
    for it = 1:numbofparts
        % Sample vn
        
        if ( i == 1 )
            prob = 1 / m * ones(1, m);
            newtheta = mu + sqrt((1 - a) * sigmasq) * randn(1, m);
            logprob = - 0.5 * (data(i) - newtheta).^2 / (a * sigmasq);
            fprob = cumsum(prob .* exp(logprob - max(logprob)));
            logweight(it) = log(fprob(end)) + max(logprob);
            fprob = fprob / fprob(end);
            u1 = rand;
            sstar = 1;
            while ( fprob(sstar) < u1 )
                sstar = sstar + 1;
            end
            
            theta{1, it} = newtheta(sstar);
            sumy{1, it} = data(i);
            sumysq{1, it} = data(i)^2;
            nj{1, it} = 1;
            s(it, i) = 1;
        else
            prob = [nj{1, it} M/m*ones(1, m)] / (sum(nj{1, it}) + M);
            newtheta = mu + sqrt((1 - a) * sigmasq) * randn(1, m);
            
            logprob = zeros(1, length(theta{1, it}) + m);
            logprob(1:length(theta{1, it})) = - 0.5 * (data(i) - theta{1, it}).^2 / (a * sigmasq);
            logprob((length(theta{1, it}) + 1):(length(theta{1, it}) + m)) = - 0.5 * (data(i) - newtheta).^2 / (a * sigmasq);
            fprob = cumsum(prob .* exp(logprob - max(logprob)));
            logweight(it) = log(fprob(end)) + max(logprob);
            fprob = fprob / fprob(end);
            u1 = rand;
            sstar = 1;
            while ( fprob(sstar) < u1 )
                sstar = sstar + 1;
            end
            if ( sstar <= length(nj{1, it}) )
                s(it, i) = sstar;
                nj{1, it}(s(it, i)) = nj{1, it}(s(it, i)) + 1;
                sumy{1, it}(s(it, i)) = sumy{1, it}(s(it, i)) + data(i);
                sumysq{1, it}(s(it, i)) = sumysq{1, it}(s(it, i)) + data(i)^2;
            else
                s(it, i) = length(nj{1, it}) + 1;
                theta{1, it} = [theta{1, it} newtheta(sstar - length(nj{1, it}))];
                nj{1, it} = [nj{1, it} 1];
                sumy{1, it} = [sumy{1, it} data(i)];
                sumysq{1, it} = [sumysq{1, it} data(i)^2];
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
    theta = theta(1, partstar);
    s = s(partstar, :);        
end