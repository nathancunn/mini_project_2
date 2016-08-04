alldata = importdata('Data/GaussianTestData1.csv', ',',1);
dataForThisContext = alldata.data;
nGenes = length(dataForThisContext);
nFeatures = size(dataForThisContext, 2);
% Set a new s for the particle filter called spart
spart = zeros(numbofparts, nGenes);
logweight = zeros(1, numbofparts);
partstar = zeros(1, numbofparts);
M = 1;
a = 0.5;
N = 4;
mu = cell(1, numbofparts);
mu(1, :) = {zeros(N, nFeatures)};
sigmasq = cell(1, numbofparts);
sigmasq(1, :) = {1 + zeros(N, nFeatures)};

% Use nGenes for n
% Make these greater dimension to allow for the multiple datasets; not sure yet
sumy = cell(1, numbofparts);
sumysq = cell(1, numbofparts);
nj = cell(1, numbofparts);
sumv = zeros(numbofparts, 1);
% Some specifications for the finite mixture model
sumy(1, :) = {zeros(N, nFeatures)};
sumysq(1, :) = {zeros(N, nFeatures)};
nj(1, :) = {zeros(N, 1)};

for i = 1:nGenes
    
    for part = 1:numbofparts
        % Sample vn - not needed?
        % sumv(part) = sumv(part) + samplev_NGG(NGGlambda, sumv(part), i, length(nj{1, part}), M, gamma);
        if ( i == 1 ) % If we're at the first data point
            sumy{1, part}(1, :) = dataForThisContext(i, :);
            nj{1, part}(1) = 1;
            logweight(part) = 0;
            spart(part, i) = 1;
        else
            % Prob uses the prob calculated previously
            % prob = [(nj{1, part}-gamma)  (M*(NGGlambda + sumv(part))^gamma)];
            prob = prob / sum(prob); % Normalising
            divisor = repmat(nj{1, part} / a + 1 / (1 - a), 1, nFeatures);
            mustar = (sumy{1, part} / a + mu{1, part} / (1 - a)) ./ divisor;
            varstar = sigmasq{1, part} ./ divisor;
            logprob = zeros(N, nFeatures);
            for g  = 1:N
                logprob(g, :) = - 0.5 * (dataForThisContext(i, :) - mustar(g, :)).^2 ./ (a * sigmasq{1, part}(g, :) + varstar(g, :)) - 0.5 * log(a * sigmasq{1, part}(g, :) + varstar(g, :));
            end
            %logprob = - 0.5 * (data(i, :) - mustar).^2 ./ (a * sigmasq + varstar) - 0.5 * log(a * sigmasq + varstar);
            logprob = transpose(log(sum(exp(logprob), 2)));
            fprob = cumsum(prob .* exp(logprob - max(logprob)));
            logweight(part) = logweight(part) + log(fprob(end)) + max(logprob);
            fprob = fprob / fprob(end);
            u1 = rand;
            sstar = 1;
            while ( fprob(sstar) < u1 )
                sstar = sstar + 1;
            end
            
            spart(part, i) = sstar;
            nj{1, part}(spart(part, i)) = nj{1, part}(spart(part, i)) + 1;
            sumy{1, part}(spart(part, i), :) = sumy{1, part}(spart(part, i), :) + dataForThisContext(i, :);
        end
    end
    
    
    ESS = sum(exp(logweight - max(logweight)))^2 / sum(exp(logweight - max(logweight)).^2);
    % This to be done at all steps?
    if ( ESS < 0.5 * numbofparts )
        fprob = cumsum(exp(logweight - max(logweight)));
        fprob = fprob / fprob(end);
        
        u1 = rand / numbofparts;
        m = 1;
        it = 1;
        while ( m <= numbofparts )
            while ( (u1 < fprob(m)) && (m <= numbofparts) )
                partstar(it) = m;
                u1 = u1 + 1 / numbofparts;
                it = it + 1;
            end
            m = m + 1;
        end
        
        sumy = sumy(1, partstar);
        nj = nj(1, partstar);
        sumv = sumv(partstar);
        spart = spart(partstar, :);
        
        logweight = zeros(1, numbofparts);
        
        for it = 1:numbofparts
            for m = 1:(i-1)
                nj{1, it}(spart(it, m)) = nj{1, it}(spart(it, m)) - 1;
                sumy{1, it}(spart(it, m), :) = sumy{1, it}(spart(it, m), :) - dataForThisContext(j, :);
                
                % This line isn't needed
                if ( nj{1, it}(spart(it, j)) == 0 )
                    nj{1, it}(spart(it, j)) = [];
                    sumy{1, it}(spart(it, j)) = [];
                    spart(it, :) = spart(it, :) - (spart(it, :) > spart(it, j));
                end
                
                
                divisor = repmat(nj{1, part} / a + 1 / (1 - a), 1, nFeatures);
                mustar = (sumy{1, part} / a + mu{1, part} / (1 - a)) ./ divisor;
                varstar = sigmasq{1, part} ./ divisor;
                logprob = zeros(N, nFeatures);
                for g  = 1:N
                    logprob(g, :) = - 0.5 * (dataForThisContext(i, :) - mustar(g, :)).^2 ./ (a * sigmasq{1, part}(g, :) + varstar(g, :)) - 0.5 * log(a * sigmasq{1, part}(g, :) + varstar(g, :));
                end
                %logprob = - 0.5 * (data(i, :) - mustar).^2 ./ (a * sigmasq + varstar) - 0.5 * log(a * sigmasq + varstar);
                logprob = transpose(log(sum(exp(logprob), 2)));
                prob = prob .* exp(- 0.5 * (dataForThisContext(j) - mustar).^2 ./ (a * sigmasq + varstar) - 0.5 * log(a * sigmasq + varstar));
                
                fprob = cumsum(prob);
                fprob = fprob / fprob(end);
                u1 = rand;
                sstar = 1;
                while ( fprob(sstar) < u1 )
                    sstar = sstar + 1;
                end
                spart(part, i) = sstar;
                nj{1, part}(spart(part, i)) = nj{1, part}(spart(part, i)) + 1;
                sumy{1, part}(spart(part, i), :) = sumy{1, part}(spart(part, i), :) + dataForThisContext(i, :);
            end
        end
    elseif (i == nGenes )
        fprob = cumsum(exp(logweight - max(logweight)));
        fprob = fprob / fprob(end);
        u1 = rand / numbofparts;
        m = 1;
        it = 1;
        while ( m <= numbofparts )
            while ( (u1 < fprob(m)) && (m <= numbofparts) )
                partstar(it) = m;
                u1 = u1 + 1 / numbofparts;
                it = it + 1;
            end
            m = m + 1;
        end
        
        sumy = sumy(1, partstar);
        nj = nj(1, partstar);
        sumv = sumv(partstar);
        spart = spart(partstar, :);
    end
    
end