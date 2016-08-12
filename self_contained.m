outputDir = '';
numberOfComponents = 5;
numbofits = 100;
fileNames = {'Data/synth_data/testData1.csv', 'Data/synth_data/testData2.csv'};
%fileNames = {'Data/GaussianTestData1.csv', 'Data/GaussianTestData2.csv'};
%fileNames = {'Data/synth_data/GaussianTestData1.csv', 'Data/synth_data/GaussianTestData2.csv'};
dataTypes = {'Gaussian', 'Gaussian'};
hyperSamplingFrequency = 1;
samplingFrequency = 10;
uniqueIdentifier = 2;
initialise = true;

global K N twoPowKminus1 powersOfTwo contextInvolvementIndexMatrix
global twoPowNumberOfPhis binaryMatrixNumberOfPhis phiIndexMatrix
global uniqueCoefficientIndexMatrix nUniqueCoefficients fHandles
global binaryMatrixKWithoutFirstRow  numberOfPhis nGenes oneMatrix
global doNotPertainToContexti J contextColumnProductIndexMatrix
global allLabelsMatrix finalIndexMatrix twoPowKminusKminus1 columnProductIndexMatrix
global timeCourseSwitches gaussianSwitches poissonSwitches nbSwitches
%%----------------------------------------------------------------------
%% SEED THE PSEUDO-RANDOM NUMBER GENERATOR -----------------------------
%%----------------------------------------------------------------------
seed1 = sum(100*clock) + sum(100*double(uniqueIdentifier));
seed2 = sum(100*clock) + sum(100*double(uniqueIdentifier));
randn('seed', seed1);%%clock-seed the random number generator (with a chain-depenedent offset)
rand('seed', seed2);%%clock-seed the random number generator (with a chain-depenedent offset)
%%----------------------------------------------------------------------
%% CONSTRUCT A RUN NAME ------------------------------------------------
%%----------------------------------------------------------------------
remains  = strtok(fileNames, '.');
while any(strcmp(remains, '')==0)
    [shortFileNames, remains] = strtok(remains, '/');
end
runName = shortFileNames{1};
for i=2:length(shortFileNames)
    runName = [runName, '_', shortFileNames{i}];
end
%%----------------------------------------------------------------------
%% OPTIONAL FEATURE SELECTION ------------------------------------------
%%----------------------------------------------------------------------
%%if not given the switches as an input, we define the feature
%selection to be switched off.
%%Note that feature selection is only acted upon by the Multinomial
%and Gaussian data types
if exist('featureSelectionSwitches', 'var')==0
    featureSelectionSwitches = zeros(1, length(fileNames));
end
featureSelectionSwitches = featureSelectionSwitches(:);%%standardise the format
%%----------------------------------------------------------------------
%% FIND USEFUL VALUES --------------------------------------------------
%%----------------------------------------------------------------------
numberOfDataSets  = length(fileNames);
outputPath        = [outputDir, runName, '_', num2str(uniqueIdentifier)];
mcmcFile          = [outputPath, '_mcmcSamples.csv'];
biomarkerFile     = [outputPath, '_featureParameters.csv'];
splitMergeCounter = zeros(1,2);
%%----------------------------------------------------------------------
%% INITIALISE PARAMETERS -----------------------------------------------
%%----------------------------------------------------------------------
K                            = numberOfDataSets;
N                            = numberOfComponents;
allLabelsMatrix              = (1:N)'*ones(1,K);
numberOfPhis                 = K*(K-1)/2;
twoPowK                      = 2^K;
twoPowKminus1                = twoPowK-1;
twoPowKminusKminus1          = twoPowKminus1 - K;
oneMatrix                    = ones(2^(K*(K-1)/2),1);
binaryMatrixK                = dec2bin(0:2^K-1) - '0';
binaryMatrixKWithoutFirstRow = binaryMatrixK(2:end,:);
powersOfTwo                  = 2.^(0:K-1);
twoPowNumberOfPhis           = 2^numberOfPhis;
binaryMatrixNumberOfPhis     = logical(dec2bin(0:twoPowNumberOfPhis-1) - '0');
% We will keep a track of the hyperparameter acceptances
nHyperProposalsMatrix   = zeros(K,3);
nHyperAcceptancesMatrix = zeros(K,3);
massParam               = 2+zeros(1,K);
gammaMatrix             = zeros(N, K);
% Initialise the gammaMatrix by sampling from the prior
for i = 1:K
    gammaMatrix(:, i) = gamrnd(massParam(i) / N, 1,N,1) + realmin;
end
%%DEPENDENCE PARAMETERS
%Sample from a Gamma(1, 0.2) prior
phiVector = gamrnd(1,1/.2, 1, numberOfPhis);

timeCourseSwitches    = false(K,1);
multinomialSwitches   = false(K,1);
bagOfWordsSwitches    = false(K,1);
gaussianSwitches      = false(K,1);% PDWK 20120120 - added gaussianSwitches
poissonSwitches       = false(K,1);
nbSwitches            = false(K,1);
fHandles              = cell(K,1);
clusterContainer(1,K) = struct('clusterStruct',[], ...
    'data', [],...
    'geneNames', [],...
    'featureNames', [], ...
    'nGenes', [], ...
    'nFeatures', [], ...
    'maxNumberOfComponents', [], ...
    'featureSelectionSwitch', []');
% Read in the data for the first context, in order to determine the number
% of genes:
dataType = dataTypes{1};
fileName = fileNames{1};
% set function handles:
switch dataType
    case 'TimeCourseEqualSpacing'
        fHandles{1}           = @TimeCourseEqualSpacing;
        timeCourseSwitches(1) = true;
    case 'TimeCourseUnequalSpacing'
        fHandles{1}           = @TimeCourse;
        timeCourseSwitches(1) = true;
    case 'Multinomial'
        fHandles{1}            = @Multinomial;
        multinomialSwitches(1) = true;
    case 'BagOfWords'
        fHandles{1}           = @BagOfWords;
        bagOfWordsSwitches(1) = true;
    case 'Gaussian' % PDWK 20120120 - added Gaussian data type
        fHandles{1}         = @Gaussian;
        gaussianSwitches(1) = true;
    case 'Poisson'
        fHandles{1}        = @Poisson;
        poissonSwitches(1) = true;
    case 'NegativeBinomial'
        fHandles{1}         = @NegativeBinomial;
        nbSwitches(1)       = true;
end
% Read in data (genenames, featurenames, and the data)
allData = importdata(fileName, ',',1);
data    = allData.data;
if(strcmp(dataType, 'Multinomial'))
    %We require the data to be numbers 1,2, ..., nLevels
    dataLevels = unique(data)';
    nLevels    = length(dataLevels);
    if(~isequal(dataLevels, 1:nLevels))
        fprintf('relabeling multinomial data:\n')
        disp(dataLevels)
        fprintf('to\n')
        disp(requiredLevels)
        
        %Force the required format
        old = data;
        for i = 1:nLevels
            currentLevel              = dataLevels(i);
            data(old == currentLevel) = i;
        end
        clear old
    end
end
%%IF GAUSSIAN DATA, PERFORM A SIMPLE OVERALL NORMALISATION
%%this ensures assumptions in our data model are reasonable
%%use robust values to do the normalisation, as data may contain
%high value outliers
if strcmp(dataType, 'Gaussian')
    sigmaValue = 0.5 * (quantile(data(:), 0.5) - quantile(data(:), 0.05));
    data       = data - median(data(:));
    data       = data / sigmaValue;
end
%%FIND/STORE VALUES FOR THE FIRST DATA TYPE
if isfield(allData, 'featureNames') && isfield(allData, 'geneNames')
    featureNames                           = allData.featureNames;
    geneNames                              = allData.geneNames;
else
    featureNames                           = allData.textdata(1,2:end);
    geneNames                              = allData.textdata(2:end,1);
end
nGenes                                     = length(geneNames);
nFeatures                                  = length(featureNames);
clusterContainer(1).data                   = data;
clusterContainer(1).geneNames              = geneNames;
clusterContainer(1).featureNames           = featureNames;
clusterContainer(1).nGenes                 = nGenes;
clusterContainer(1).nFeatures              = nFeatures;
clusterContainer(1).maxNumberOfComponents  = N;
clusterContainer(1).featureSelectionSwitch = featureSelectionSwitches(1);
clusterContainer(1).featureParameters      = ones(1, nFeatures);
%%OPTION TO RANDOMLY INITIALISE THE FEATURE PARAMETERS
if featureSelectionSwitches(1)
    clusterContainer(1).featureParameters = randi([0 1], 1, nFeatures);
end
% Initialise the clusters for the first component:
[structOfClusters clusterIDs]     = feval(fHandles{1}, clusterContainer(1), 'init');
clusterContainer(1).clusterStruct = structOfClusters;

s      = zeros(nGenes, K);
s(:,1) = clusterIDs;


% Now read in the data for the remaining contexts:
for k = 2:K
    dataType = dataTypes{k};
    fileName = fileNames{k};
    % set function handles:
    switch dataType
        case 'TimeCourseEqualSpacing'
            fHandles{k}           = @TimeCourseEqualSpacing;
            timeCourseSwitches(k) = true;
        case 'TimeCourseUnequalSpacing'
            fHandles{k}           = @TimeCourse;
            timeCourseSwitches(k) = true;
        case 'Multinomial'
            fHandles{k}            = @Multinomial;
            multinomialSwitches(k) = true;
        case 'BagOfWords'
            fHandles{k}           = @BagOfWords;
            bagOfWordsSwitches(k) = true;
        case 'Gaussian' % PDWK 20120120 - added Gaussian data type
            fHandles{k}         = @Gaussian;
            gaussianSwitches(k) = true;
        case 'Poisson'
            fHandles{k}        = @Poisson;
            poissonSwitches(k) = true;
        case 'NegativeBinomial'
            fHandles{k}         = @NegativeBinomial;
            nbSwitches(k)       = true;
    end
    % Read in data (genenames, featurenames, and the data)
    allData = importdata(fileName, ',',1);
    data    = allData.data;
    if(strcmp(dataType, 'Multinomial'))
        %We require the data to be numbers 1,2, ..., nLevels
        dataLevels = unique(data)';
        nLevels    = length(dataLevels);
        if(~isequal(dataLevels, 1:nLevels))
            fprintf('relabeling multinomial data:\n')
            disp(dataLevels)
            fprintf('to\n')
            disp(requiredLevels)
            
            %Force the required format
            old = data;
            for i = 1:nLevels
                currentLevel              = dataLevels(i);
                data(old == currentLevel) = i;
            end
            clear old
        end
    end
    %%IF GAUSSIAN DATA, PERFORM A SIMPLE OVERALL NORMALISATION
    %%this ensures assumptions in our data model are reasonable
    %%use robust values to do the normalisation, as data may contain
    %high value outliers
    if strcmp(dataType, 'Gaussian')
        sigmaValue = 0.5 * (quantile(data(:), 0.5) - quantile(data(:), 0.05));
        data       = data - median(data(:));
        data       = data / sigmaValue;
    end
    %%FIND/STORE VALUES FOR THE K-th DATA TYPE
    if isfield(allData, 'featureNames') && isfield(allData, 'geneNames')
        featureNames                           = allData.featureNames;
        geneNames                              = allData.geneNames;
    else
        featureNames                           = allData.textdata(1,2:end);
        geneNames                              = allData.textdata(2:end,1);
    end
    nGenes                                     = length(geneNames);
    nFeatures                                  = length(featureNames);
    clusterContainer(k).data                   = data;
    clusterContainer(k).geneNames              = geneNames;
    clusterContainer(k).featureNames           = featureNames;
    clusterContainer(k).nGenes                 = nGenes;
    clusterContainer(k).nFeatures              = nFeatures;
    clusterContainer(k).maxNumberOfComponents  = N;
    clusterContainer(k).featureSelectionSwitch = featureSelectionSwitches(k);
    clusterContainer(k).featureParameters      = ones(1, nFeatures);
    %%OPTION TO RANDOMLY INITIALISE THE FEATURE PARAMETERS
    if featureSelectionSwitches(k)
        clusterContainer(k).featureParameters = randi([0 1], 1, nFeatures);
    end
    %%INITIALISE THE CLUSTERS
    [structOfClusters clusterIDs]     = feval(fHandles{k}, clusterContainer(k), 'init');
    clusterContainer(k).clusterStruct = structOfClusters;
    s(:,k)                            = clusterIDs;
end

for k = 1:K
    proposedClustersForThisContext = clusterContainer(k).clusterStruct;
    dataForThisContext = clusterContainer(k).data;
    for ind = 1:nGenes
        if(s(ind, k) <= numberOfComponents)
        proposedClustersForThisContext = AddRemoveItem('removeGene', proposedClustersForThisContext, s(ind, k), ind, dataForThisContext, 1);
        end
    end
    clusterContainer(k).clusterStruct = proposedClustersForThisContext;
end
numbofparts = 10;
proposedClusterContainer = repelem(clusterContainer, numbofparts);

s = zeros(nGenes, K, numbofparts);
logweight = zeros(K, numbofparts);

prior = transpose(gammaMatrix);




for i = 1:nGenes
    %disp(['i = ' num2str(i)])
    if(i == 1)
        fprintf ('i = %d,', i)
    elseif(mod(i, 20) ==0)
        fprintf (' %d,\r   ', i)
    else
        fprintf (' %d,', i)
    end

            
    for m = 1:numbofparts
        for k = 1:K;
            dataForThisContext = proposedClusterContainer(numbofparts * (k - 1) + m).data;
            proposedClustersForThisContext = proposedClusterContainer(numbofparts * (k - 1) + m).clusterStruct;
            if(i == 1)
                proposedClustersForThisContext = AddRemoveItem('addGene', proposedClustersForThisContext, 1, i, dataForThisContext, 1);
                proposedClusterContainer(numbofparts * (k - 1) + m).clusterStruct = proposedClustersForThisContext;
                sstar = 1;
            else
                logprob = zeros(1, numberOfComponents);
                oldLogMarginalLikelihoods         = [proposedClustersForThisContext.logMarginalLikelihood];
                for ind = 1:numberOfComponents
                    % prop = AddRemoveItem('addGene', proposedClustersForThisContext, ind, i, dataForThisContext, 1);
                    %if (proposedClustersForThisContext(ind).nGenes > 1)
                    %    logprob(1, ind) = mvnpdf(dataForThisContext(i, :), proposedClustersForThisContext(ind).empiricalMeans, cov(dataForThisContext(find(proposedClustersForThisContext(ind).logicalGeneIDs == 1), :)));
                    %else
                        logprob(1, ind) = mvnpdf(dataForThisContext(i, :), proposedClustersForThisContext(ind).empiricalMeans);
                    %end
                    %proposedClustersForThisContext = AddRemoveItem('addGene', proposedClustersForThisContext, ind, i, dataForThisContext, 1);
                    % newLogMarginalLikelihoods         = [prop.logMarginalLikelihood];
                    %marginalLikelihoodRatio           = newLogMarginalLikelihoods - oldLogMarginalLikelihoods;
                    %logprob(1, ind) = sum(newLogMarginalLikelihoods) - oldLogMarginalLikelihoods;
                    % logprob(1, ind) = newLogMarginalLikelihoods(ind); % - oldLogMarginalLikelihoods(ind);
                end
                %newLogMarginalLikelihoods         = [proposedClustersForThisContext.logMarginalLikelihood];
                % marginalLikelihoodRatio           = newLogMarginalLikelihoods - oldLogMarginalLikelihoods;
                %marginalLikelihoodRatio = marginalLikelihoodRatio - max(marginalLikelihoodRatio);
                % logprob = marginalLikelihoodRatio;
                %logprob = newLogMarginalLikelihoods;
                %logprob(1,ind) = (sum(([prop.logMarginalLikelihood])));
                prob = prior(k, :)/sum(prior(k, :));
                % prob = ones(1, numberOfComponents);
                fprob = cumsum(prob .* exp(logprob - max(logprob)));
                fprob = fprob/fprob(end);
                u1 = rand;
                sstar = 1;
                while ( fprob(sstar) < u1 )
                    sstar = sstar + 1;
                end
                logweight(k, m) = logweight(k, m) + log(fprob(end)) + max(logprob);
                proposedClusterContainer(numbofparts * (k - 1) + m).clusterStruct = AddRemoveItem('addGene', proposedClusterContainer(numbofparts * (k - 1) + m).clusterStruct, sstar, i, dataForThisContext, 1);
                
            end
            s(i, k, m) = sstar;
            
        end
    end
end


ties = find(sum(logweight) == max(sum(logweight)));
partstar = ties(randi(length(ties)));
s = s(:, :, partstar)
tabulate(s(:, 1))
for k = 1:K
    clusterContainer(k) = proposedClusterContainer(numbofparts * (k - 1) + partstar)
end

prob     = gammaMatrix(:, k);%%conditional prior
%%We need to find the labels of gene i in the other contexts:
%%We consider all possibilities for the label of gene i
%%in context j:
labelsAcrossAllContexts            = s(i, :, m);
labelsAcrossAllContextsMatrix      = ones(N,1)*labelsAcrossAllContexts;
labelsAcrossAllContextsMatrix(:,k) = 1:N;
%%We then need to upweight the probabilities according
%%to the labels in the other contexts.
labelAgreementMatrix = (allLabelsMatrix == labelsAcrossAllContextsMatrix);
myPhiMatrix          = ones(size(labelAgreementMatrix,1),1)*phiVector;
binInd               = labelAgreementMatrix*(2.^(size(labelAgreementMatrix,2)-1:-1:0))';  %bin2dec(num2str(labelAgreementMatrix))
finalIndexMatrixRows = finalIndexMatrix(binInd,:);
finalIndexMatrixRows(:,notPertinentInThisContext) = false;

upWeighting          = finalIndexMatrixRows.*myPhiMatrix;
upWeighting          = prod(1+ upWeighting(:,phiIndicesForConditional),2);
prob                 = prob.*upWeighting;  %%This takes care of multiplying by 1+phi
prob = transpose(prob);
