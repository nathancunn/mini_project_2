function [s] = MDI_CPF(s_old, data, mu, sigmasq, a, M, numbofparts)

% Read in global variables
global K nGenes N phiIndexMatrix finalIndexMatrix fHandles ...
    doNotPertainToContexti allLabelsMatrix timeCourseSwitches ...
    gaussianSwitches poissonSwitches nbSwitches massParam nComponents
allComponents = 1:N;



% Specify some values for the particle filter
sumy = cell(1, K, numbofparts);
nj = cell(1, K, numbofparts);
s = zeros(n, K, numbofparts);

logweight = zeros(1, numbofparts);
indx = zeros(1, numbofparts);

% Initialise for the first observation
for k = 1:K
    for m = 1:numbofparts
        sumy{1, k, m} = dataForThisContext(i, :);
        nj{1, k, m} = 1;
        logweight(m) = 0;
        s(1, k, m) = 1;
    end
end

% Run for all other observations
for i = 2:nGenes
    weight = exp(logweight - max(logweight));
    weight = weight / sum(weight);
    
    % Update
    fprob = cumsum(weight);
    fprob = fprob / fprob(end);
    
    for part = 1:(numbofparts - 1) % Hold last trajectory constant
        u = rand;
        indx(part) = 1;
        while( u > fprob(indx(part)))
            indx(part) = indx(part) + 1;
        end
    end
    indx(numbofparts = numbofparts;
    s(1:(i - 1), :, :) = s(1:(i - 1), :, indx);
    nj = nj(:, :, indx);
    sumy = sumy(:, :, indx);
    logweight = zeros(1, numbofparts);

    for m = 1:(numbofparts - 1)
        for k = 1:K
            dataForThisContext = clusterContainer(k).data;
            prob = [nj{1, k, m} massParam] / (sum(nj{1, k, m}) + M);
            mustar = [(sumy{1, m, k} / a + mu / (1 - a)) ./ (nj{1, m, k} / a + 1 / (1 - a)) mu];
            varstar = [sigmasq ./ (nj{1, it} / a + 1 / (1 - a)) sigmasq * (1 - a)];
            logprob = -  0.5 * (dataForThisContext(i, :) - mustar).^2 ./ (a * sigmasq + varstar) - 0.5 * log(a * sigmasq + varstar);
            fprob = cumsum(prob .* exp(logprob - max(logprob)));
            logweight(it) = logweight(it) + log(fprob(end)) + max(logprob);
            fprob = fprob / fprob(end);
            u1 = rand;
            s(i, k, m) = 1;
            while(fprob(s(i, k, m)) < u1)
                s(i, k, m) = s(i, k, m) + 1
            end
        end
    end
    
    % Re-specify the last trajectory
    for k = 1:K
        dataForThisContext = clusterContainer(k).data;
        s(i,k, numbofparts) = s_old(i, k);
        prob = [nj{1, k, numbofparts} M] / (sum(nj{1, k, numbofparts}) +  M); 
        mustar = [(sumy{1, k, numbofparts} / a + mu / (1 - a)) ./ (nj{1, k, numbofparts} / a + 1 / (1 - a)) mu];
        varstar = [sigmasq ./ (nj{1, numbofparts} / a + 1 / (1- a)) sigmasq * (1 - a)];
        logprob =  -0.5 * (dataForThisContext(i, :) - mustar) .^2 ./ (a * sigmasq + varstar) - 0.5 * log(a * sigmasq + varstar);
        fprob = cumsum(prob .* exp(logprob - max(logprob)));
        logweight(numbofparts) = logweight(numbofparts) + log(fprob(end)) + max(logprob);
    end
    
    for m = 1:numbofparts
        for k = 1:K
            dataForThisContext = clusterContainer(k).data;
            if(s(i, k, m) > length(nj{i, k, m})) % If we're starting a new cluster
                nj{1, k, m}(s(i, k, m) = 1;
                sumy{1, k, m}(s(i, k, m)) = dataForThisContext(i, :);
            else
                nj{1, k, m}(s(i, k, m)) = nj{1, k, m}(s(i, k, m)) + 1;
                sumy{1, k, m}(s(i, k, m)) = sumy{1, k, m}(s(i, k, m)) + dataForThisContext(i, :);
            end
        end
    end
end

fprob = cumsum(exp(logweight - max(logweight)));
fprob = fprob / fprob(end);

u = rand;
counter = 1;
while(u > fprob(counter))
    counter = counter + 1;
end
s = s(:, :, counter);

                
    
    




    
    

    






