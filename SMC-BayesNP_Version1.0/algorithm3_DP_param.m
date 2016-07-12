function [s, a] = algorithm3_DP_param(data, mu, sigmasq, M, numbofparts)

n = length(data);

aprior = [1 1];

sumy = cell(1, numbofparts);
sumysq = cell(1, numbofparts);
nj = cell(1, numbofparts);
sumv = zeros(numbofparts, 1);
s = zeros(numbofparts, n);
a = betarnd(aprior(1), aprior(2), numbofparts, 1);

aaccept = 0;
acount = 0;
logasd = log(0.01);

logweight = zeros(1, numbofparts);
partstar = zeros(1, numbofparts);

for i = 1:n
    disp(['i = ' num2str(i)]);
    
    for it = 1:numbofparts
        % Sample vn
        lambda = gamrnd(M + i - 1, 1 / (1 + sumv(it)));
        sumv(it) = sumv(it) - log(rand) / lambda;
        
        if ( i == 1 )
            sumy{1, it} = data(i);
            sumysq{1, it} = data(i)^2;
            nj{1, it} = 1;
            logweight(it) = 0;
            s(it, i) = 1;
        else
            prob = [nj{1, it} M] / (sum(nj{1, it}) + M);
            mustar = [(sumy{1, it} / a(it) + mu / (1 - a(it))) ./ (nj{1, it} / a(it) + 1 / (1 - a(it))) mu];
            varstar = [sigmasq ./ (nj{1, it} / a(it) + 1 / (1 - a(it))) sigmasq*(1-a(it))];
            logprob = - 0.5 * (data(i) - mustar).^2 ./ (a(it) * sigmasq + varstar) - 0.5 * log(a(it) * sigmasq + varstar);
            fprob = cumsum(prob .* exp(logprob - max(logprob)));
            logweight(it) = logweight(it) + log(fprob(end)) + max(logprob);
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
        sumv = sumv(partstar);
        s = s(partstar, :);
        a = a(partstar);
        
        logweight = zeros(1, numbofparts);
        
        for it = 1:numbofparts
            for j = 1:(i-1)
                nj{1, it}(s(it, j)) = nj{1, it}(s(it, j)) - 1;
                sumy{1, it}(s(it, j)) = sumy{1, it}(s(it, j)) - data(j);
                sumysq{1, it}(s(it, j)) = sumysq{1, it}(s(it, j)) - data(j)^2;
                
                if ( nj{1, it}(s(it, j)) == 0 )
                    nj{1, it}(s(it, j)) = [];
                    sumy{1, it}(s(it, j)) = [];
                    sumysq{1, it}(s(it, j)) = [];
                    s(it, :) = s(it, :) - (s(it, :) > s(it, j));
                end
                
                prob = [nj{1, it} M] / (sum(nj{1, it}) + M);
                mustar = [(sumy{1, it} / a(it) + mu / (1 - a(it))) ./ (nj{1, it} / a(it) + 1 / (1 - a(it))) mu];
                varstar = [sigmasq ./ (nj{1, it} / a(it) + 1 / (1 - a(it))) sigmasq*(1-a(it))];
                prob = prob .* exp(- 0.5 * (data(j) - mustar).^2 ./ (a(it) * sigmasq + varstar) - 0.5 * log(a(it) * sigmasq + varstar));
                fprob = cumsum(prob);
                fprob = fprob / fprob(end);
                u1 = rand;
                s(it, j) = 1;
                while ( fprob(s(it, j)) < u1 )
                    s(it, j) = s(it, j) + 1;
                end
                if ( s(it, j) == length(nj{1, it}) + 1 )
                    nj{1, it} = [nj{1, it} 1];
                    sumy{1, it} = [sumy{1, it} data(j)];
                    sumysq{1, it} = [sumysq{1, it} data(j)^2];
                else
                    nj{1, it}(s(it, j)) = nj{1, it}(s(it, j)) + 1;
                    sumy{1, it}(s(it, j)) = sumy{1, it}(s(it, j)) + data(j);
                    sumysq{1, it}(s(it, j)) = sumysq{1, it}(s(it, j)) + data(j)^2;
                end
            end
        end
        
        
        for it = 1:numbofparts
            
            % Update a
            mustar = (sumy{1, it} / a(it) + mu / (1 - a(it))) ./ (nj{1, it} / a(it) + 1 / (1 - a(it)));
            varstar = sigmasq ./ (nj{1, it} / a(it) + 1 / (1 - a(it)));
            theta = mustar + sqrt(varstar) .* randn(1, length(mustar));
            
            alphastar = aprior(1) - 0.5 * i;
            betastar = aprior(2) - 0.5 * length(nj{1, it});
            a1 = sum(sumysq{1, it} - 2 * sumy{1, it} .* theta + nj{1, it} .* theta.^2) / sigmasq;
            b1 = sum((theta - mu).^2) / sigmasq;
            
            %         if ( alphastar + betastar < 0 )
            %             a(it) = Updatea(alphastar, betastar, a1, b1);
            %         else
            trans = log(a(it)) - log(1 - a(it));
            newtrans = trans + exp(logasd) * randn;
            newa = exp(newtrans) / (1 + exp(newtrans));
            
            loglike = (alphastar - 1) * log(a(it)) + (betastar - 1) * log(1 - a(it));
            loglike = loglike - 0.5 * b1 / (1 - a(it)) - 0.5 * a1 / a(it);
            
            newloglike = (alphastar - 1) * log(newa) + (betastar - 1) * log(1 - newa);
            newloglike = newloglike - 0.5 * b1 / (1 - newa) - 0.5 * a1 / newa;
            
            logaccept = newloglike - loglike + log(1 / a(it) + 1 / (1 - a(it))) - log(1 / newa + 1 / (1 - newa));
            
            accept = 1;
            if (logaccept < 0 )
                accept = exp(logaccept);
            end
            
            aaccept = aaccept + accept;
            acount = acount + 1;
            
            
            if ( rand < accept )
                a(it) = newa;
            end
            
            logasd = logasd + 1 / it^0.55 * (accept - 0.3);
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
        sumv = sumv(partstar);
        s = s(partstar, :);
        a = a(partstar);        
    end
end

