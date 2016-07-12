function [s] = algorithm4_NGG(data, mu, sigmasq, a, gamma, M, m, numbofparts)

n = length(data);

NGGlambda = 1;

sumy = cell(1, numbofparts);
sumysq = cell(1, numbofparts);
nj = cell(1, numbofparts);
sumv = zeros(numbofparts, 1);
s = zeros(numbofparts, n);
theta = cell(1, numbofparts);

logweight = zeros(1, numbofparts);
partstar = zeros(1, numbofparts);

for i = 1:n
    disp(['i = ' num2str(i)]);
    
    for it = 1:numbofparts
        % Sample vn
        sumv(it) = sumv(it) + samplev_NGG(NGGlambda, sumv(it), i, length(nj{1, it}), M, gamma);
        
        if ( i == 1 )
            prob = 1 / m * ones(1, m);
            newtheta = mu + sqrt((1 - a) * sigmasq) * randn(1, m);
            logprob = - 0.5 * (data(i) - newtheta).^2 / (a * sigmasq) - 0.5 * log(a * sigmasq);
            fprob = cumsum(prob .* exp(logprob - max(logprob)));
            logweight(it) = log(fprob(end)) + max(logprob);
            fprob = fprob / fprob(end);
            u1 = rand;
            sstar = 1;
            while ( fprob(sstar) < u1 )
                sstar = sstar + 1;
            end
            
            logweight(it) = logweight(it) + log(fprob(end)) + max(logprob);
            theta{1, it} = newtheta(sstar);
            sumy{1, it} = data(i);
            sumysq{1, it} = data(i)^2;
            nj{1, it} = 1;
            s(it, i) = 1;
        else
            prob = [(nj{1, it}-gamma)  (M*(NGGlambda + sumv(it))^gamma)/m*ones(1, m)];
            prob = prob / sum(prob);
            newtheta = mu + sqrt((1 - a) * sigmasq) * randn(1, m);
            
            logprob = zeros(1, length(theta{1, it}) + m);
            logprob(1:length(theta{1, it})) = - 0.5 * (data(i) - theta{1, it}).^2 / (a * sigmasq) - 0.5 * log(a * sigmasq);
            logprob((length(theta{1, it}) + 1):(length(theta{1, it}) + m)) = - 0.5 * (data(i) - newtheta).^2 / (a * sigmasq) - 0.5 * log(a * sigmasq);
            fprob = cumsum(prob .* exp(logprob - max(logprob)));
            logweight(it) = logweight(it) + log(fprob(end)) + max(logprob);
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
        
        if ( length(theta{1, it}) ~= length(unique(s(it, 1:i))) )
            theta{1, it}
            s(it, 1:i)
            stop;
        end
        if ( length(theta{1, it}) ~= length(nj{1, it}) )
            stop;
        end
    end
    ESS = sum(exp(logweight - max(logweight)))^2 / sum(exp(logweight - max(logweight)).^2);
    
    if ( ESS < 0.5 * numbofparts )
        
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
        sumv = sumv(partstar);
        s = s(partstar, :);
        
        logweight = zeros(1, numbofparts);
        
        
        %        disp('update')
        for it = 1:numbofparts
            for j = 1:(i-1)
                nj{1, it}(s(it, j)) = nj{1, it}(s(it, j)) - 1;
                sumy{1, it}(s(it, j)) = sumy{1, it}(s(it, j)) - data(j);
                sumysq{1, it}(s(it, j)) = sumysq{1, it}(s(it, j)) - data(j)^2;
                
                newtheta = mu + sqrt((1 - a) * sigmasq) * randn(1, m);
                if ( nj{1, it}(s(it, j)) == 0 )
                    newtheta(1) = theta{1, it}(s(it, j));
                    nj{1, it}(s(it, j)) = [];
                    sumy{1, it}(s(it, j)) = [];
                    sumysq{1, it}(s(it, j)) = [];
                    theta{1, it}(s(it, j)) = [];
                    s(it, :) = s(it, :) - (s(it, :) > s(it, j));
                end
                
                prob = [(nj{1, it}-gamma)  (M*(NGGlambda + sumv(it))^gamma)/m*ones(1, m)];
                
                logprob = zeros(1, length(theta{1, it}) + m);
                logprob(1:length(theta{1, it})) = - 0.5 * (data(j) - theta{1, it}).^2 / (a * sigmasq);
                logprob((length(theta{1, it}) + 1):(length(theta{1, it}) + m)) = - 0.5 * (data(j) - newtheta).^2 / (a * sigmasq);
                fprob = cumsum(prob .* exp(logprob - max(logprob)));
                fprob = fprob / fprob(end);
                u1 = rand;
                sstar = 1;
                while ( fprob(sstar) < u1 )
                    sstar = sstar + 1;
                end
                if ( sstar <= length(nj{1, it}) )
                    s(it, j) = sstar;
                    nj{1, it}(s(it, j)) = nj{1, it}(s(it, j)) + 1;
                    sumy{1, it}(s(it, j)) = sumy{1, it}(s(it, j)) + data(j);
                    sumysq{1, it}(s(it, j)) = sumysq{1, it}(s(it, j)) + data(j)^2;
                else
                    s(it, j) = length(nj{1, it}) + 1;
                    theta{1, it} = [theta{1, it} newtheta(sstar - length(nj{1, it}))];
                    nj{1, it} = [nj{1, it} 1];
                    sumy{1, it} = [sumy{1, it} data(j)];
                    sumysq{1, it} = [sumysq{1, it} data(j)^2];
                end
            end
            
            mustar = (sumy{1, it} / a + mu / (1 - a)) ./ (nj{1, it} / a + 1 / (1 - a));
            varstar = sigmasq ./ (nj{1, it} / a + 1 / (1 - a));
            
            theta{1, it} = mustar + sqrt(varstar) .* randn(1, length(theta{1, it}));
        end
    elseif ( i == n )
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
        sumv = sumv(partstar);
        s = s(partstar, :);
    end
end

