function [s] = algorithm3_NGG(data, mu, sigmasq, a, gamma, M, numbofparts)

n = length(data);

NGGlambda = 1;

sumy = cell(1, numbofparts);
nj = cell(1, numbofparts);
sumv = zeros(numbofparts, 1);
s = zeros(numbofparts, n);

logweight = zeros(1, numbofparts);
partstar = zeros(1, numbofparts);

for i = 1:n
    disp(['i = ' num2str(i)]);
    
    for it = 1:numbofparts
        % Sample vn
        sumv(it) = sumv(it) + samplev_NGG(NGGlambda, sumv(it), i, length(nj{1, it}), M, gamma);
        
        if ( i == 1 )
            sumy{1, it} = data(i);
            nj{1, it} = 1;
            logweight(it) = 0;
            s(it, i) = 1;
        else
            prob = [(nj{1, it}-gamma)  (M*(NGGlambda + sumv(it))^gamma)];
            prob = prob / sum(prob);

            mustar = [(sumy{1, it} / a + mu / (1 - a)) ./ (nj{1, it} / a + 1 / (1 - a)) mu];
            varstar = [sigmasq ./ (nj{1, it} / a + 1 / (1 - a)) sigmasq*(1-a)];
            logprob = - 0.5 * (data(i) - mustar).^2 ./ (a * sigmasq + varstar) - 0.5 * log(a * sigmasq + varstar);
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
            else
                s(it, i) = length(nj{1, it}) + 1;
                nj{1, it} = [nj{1, it} 1];
                sumy{1, it} = [sumy{1, it} data(i)];
            end
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
        nj = nj(1, partstar);
        sumv = sumv(partstar);
        s = s(partstar, :);
        
        logweight = zeros(1, numbofparts);
        
        for it = 1:numbofparts
            for j = 1:(i-1)
                nj{1, it}(s(it, j)) = nj{1, it}(s(it, j)) - 1;
                sumy{1, it}(s(it, j)) = sumy{1, it}(s(it, j)) - data(j);
                
                if ( nj{1, it}(s(it, j)) == 0 )
                    nj{1, it}(s(it, j)) = [];
                    sumy{1, it}(s(it, j)) = [];
                    s(it, :) = s(it, :) - (s(it, :) > s(it, j));
                end
                
                prob = [(nj{1, it}-gamma)  (M*(NGGlambda + sumv(it))^gamma)];
                mustar = [(sumy{1, it} / a + mu / (1 - a)) ./ (nj{1, it} / a + 1 / (1 - a)) mu];
                varstar = [sigmasq ./ (nj{1, it} / a + 1 / (1 - a)) sigmasq*(1-a)];
                prob = prob .* exp(- 0.5 * (data(j) - mustar).^2 ./ (a * sigmasq + varstar) - 0.5 * log(a * sigmasq + varstar));
                
                fprob = cumsum(prob);
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
                else
                    s(it, j) = length(nj{1, it}) + 1;
                    nj{1, it} = [nj{1, it} 1];
                    sumy{1, it} = [sumy{1, it} data(j)];
                end
            end
        end
    elseif (i == n )
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
        nj = nj(1, partstar);
        sumv = sumv(partstar);
        s = s(partstar, :);        
    end
end

