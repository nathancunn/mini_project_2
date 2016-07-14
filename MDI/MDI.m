%%MULTIPLE DATASET INTEGRATION (MDI)
%%----------------------------------
%%Original version by Paul Kirk
%%Also modified by Rich Savage
%%
function MDI(outputDir, numberOfComponents, numbofits, fileNames, dataTypes, hyperSamplingFrequency, samplingFrequency, uniqueIdentifier, initialise, featureSelectionSwitches)

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
% The following vector stores the component indicator variables (cluster
% IDS)
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
save([outputPath, 'SEEDS.mat'], 'seed1', 'seed2');
whichTimeCourseDataSets = find(timeCourseSwitches');
%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATE SOME INDEX MATRICES:
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%columnProductIndexMatrix -- describes the matrix columns that
%%must be multiplied together in order to find all column prods
columnProductIndexMatrix = calculateColumnProductIndexMatrix();
%%Suppose we change gammaMatrix for component n in the kth context.  Which
%%elements of columnProductMatrix must we change?  Well, we know that we
%%need only look along row n, but which columns pertain to context k?  We
%%require another index matrix to tell us:
%%contextInvolvementIndexMatrix -- for each context, tells us
%%whether or not the jth column of columnProductMatrix involves
%%context j
contextInvolvementIndexMatrix = calculateContextInvolvementIndexMatrix();
%%phiIndexMatrix -- describes the order of the phi combinations
%%that is used throughout.  e.g. (K=4) phi_12, phi_13, phi_14,
%%phi_23, phi_24, phi_34.
phiIndexMatrix = calculatePhiIndexMatrix();
%%We have a vector, phiVector, storing all of the current phis.  Which
%%elements of phiVector pertain to Context i?
doNotPertainToContexti = calculateContextPertinentPhis();
%%coefficientIndexMatrix -- for each phi, prescribes which
%%columns of sumOfColumnProducts need to be multiplied together
%%in order to find the normalizing constant
[uniqueCoefficientIndexMatrix J] = calculateCoefficientIndexMatrix();
nUniqueCoefficients              = size(uniqueCoefficientIndexMatrix,1);
%%contextColumnProductIndexMatrix -- rows indexed by contexts,
%%columns by unique coefficients.  The ij-entry tells us which
%%column of columnProductMatrix contains the contribution of
%%context i to the jth row of uniqueCoefficientIndexMatrix
contextColumnProductIndexMatrix = calculateContextColumnProductIndexMatrix();
%%finalIndexMatrix -- a simple dictionary to enable us to look
%%up the appropriate phi weighting factors for any
%%configuration of the component labels
finalIndexMatrix = calculateFinalIndexMatrix();
%% CALCULATE THE PRODUCTS OF THE COLUMNS OF THE gammaMatrix MATRIX
%%Calculate the column products for the initial gammaMatrix:
columnProductMatrix                = zeros(N,twoPowKminus1);
columnProductMatrix(:,powersOfTwo) = gammaMatrix;
% Calculate column products
for i = 1:twoPowKminusKminus1
    columnProductMatrix(:,columnProductIndexMatrix(i,3)) = columnProductMatrix(:,columnProductIndexMatrix(i,1)).*columnProductMatrix(:,columnProductIndexMatrix(i,2));
end
sumsOfColumnProducts        =  sum(columnProductMatrix);
%% CALCULATE THE NORMALIZING CONSTANT
uniquePhiCoeffVector = zeros(1,nUniqueCoefficients);
for i = 1:nUniqueCoefficients
    uniquePhiCoeffVector(i) = prod(sumsOfColumnProducts(uniqueCoefficientIndexMatrix(i,:)));
end
phiCoeffVector      = uniquePhiCoeffVector(J);
productMatrix       = phiCoeffVector'.*prod(  1+((phiVector(oneMatrix,:)-1).*binaryMatrixNumberOfPhis)  ,2);
normalizingConstant = sum(productMatrix);

if(initialise)
    %% INITIALISE THE OUTPUT FILE/S
    writeNewFile(mcmcFile, phiVector, massParam, s, geneNames)
    if any(featureSelectionSwitches)
        InitialiseBiomarkerFile(biomarkerFile, clusterContainer);
    end
else
    load([outputPath '.mat'])
end

tic()
for it = 1:numbofits
    %%----------------------------------------------------------------------
    %% PRINT SOME PROGRESS INFORMATION TO SCREEN ---------------------------
    %%----------------------------------------------------------------------
    if(mod(it, samplingFrequency) == 0)
        disp(['MCMC sample ', num2str(it)])
        toc(), tic()
        %%PLOT THE CURRENT CLUSTERING PARTITION TO SCREEN
        %    [dum, sortIndex] = sort(s(:,1));
        sortIndex = 1:size(s, 1);
        figure(1), imagesc(s(sortIndex, :)'), colorbar, xlabel('Data Items'), ylabel('Data Sets'), pause(1e-3)
        for j=1:numberOfDataSets
            figure(j+1), imagesc(clusterContainer(j).data(sortIndex, :)'), colorbar, pause(1e-3)
        end
        %%DISPLAY THE CURRENT NUMBER OF OCCUPIED CLUSTERS FOR EACH DATA TYPE
        outputString = 'Current nClusters:     ';
        for j=1:numberOfDataSets, outputString = [outputString, num2str(length(unique(s(:,j)))), ', '];, end
        outputString = outputString(1:(end-2));
        disp(outputString)
        %%OPTION TO WRITE OUT THE NUMBER OF SWITCHED-ON FEATURES
        if any(featureSelectionSwitches)
            nSwitchedOn = zeros(1, numberOfDataSets);
            for k=1:numberOfDataSets
                featureParameters = clusterContainer(k).clusterStruct.featureParameters;
                nSwitchedOn(k)    = sum(featureParameters);
            end
            disp(['nFeaturesSwitchedOn: ', num2str(nSwitchedOn)])
        end
        %%SPLIT-MERGE PERFORMANCE
        disp(['SplitMerge acceptance: ', num2str(splitMergeCounter(1)), ' of ', num2str(splitMergeCounter(2))])
    end
    %%----------------------------------------------------------------------
    %% RESAMPLE VARIOUS PARAMETERS, CLUSTER LABELS -------------------------
    %%----------------------------------------------------------------------
    vStar                      = DrawNewVStar(normalizingConstant);
    massParam                  = DrawNewMassParameter(massParam, gammaMatrix);
    [clusterContainer s]       = DrawNewItemLabel(gammaMatrix, phiVector, s, clusterContainer);
    [phiVector  productMatrix] = DrawNewPhi(phiVector, vStar, s, productMatrix);
    normalizingConstant        = sum(productMatrix);%recalculate the normalizing constant
    %%----------------------------------------------------------------------
    %% RESAMPLE THE GAMMAS -------------------------------------------------
    %%----------------------------------------------------------------------
    [gammaMatrix columnProductMatrix sumsOfColumnProducts uniquePhiCoeffVector]  = ...
        DrawNewGamma(vStar, massParam, phiVector, gammaMatrix, columnProductMatrix, sumsOfColumnProducts, uniquePhiCoeffVector, s);
    %%----------------------------------------------------------------------
    %% Allow "block updates" of the cluster labels (e.g. swap all 1s and 3s)
    %%----------------------------------------------------------------------
    [s clusterContainer gammaMatrix columnProductMatrix productMatrix sumsOfColumnProducts uniquePhiCoeffVector] = ...
        BlockLabelUpdate(s, clusterContainer, phiVector, gammaMatrix);
    %%----------------------------------------------------------------------
    %% IF REQUIRED, RESAMPLE THE FEATURE SELECTION PARAMETERS --------------
    %%----------------------------------------------------------------------
    for k=1:numberOfDataSets
        if featureSelectionSwitches(k)
            clusterContainer(k).clusterStruct = feval(fHandles{k}, clusterContainer(k).clusterStruct, 'featureSelection', s(:,k));
        end
    end
    %%----------------------------------------------------------------------
    %% SPLIT-MERGE MCMC STEP -----------------------------------------------
    %%----------------------------------------------------------------------
    %%this is a Metropolis-Hastings MCMC step, based on the
    %%split-merge algorithm of Dahl (2005):
    %%'Sequentially-allocated merge-split sampler for conjugate and
    %nonconjugate Dirichlet process mixture models'
    %%
    %%NOTE: for now, only do this for the conjugate data models
    %%the non-conjugate case is a bit more complicated
    if strcmp(dataType, 'Multinomial') | strcmp(dataType, 'Gaussian') | strcmp(dataType, 'Poisson')
        for k=1:numberOfDataSets
            for l=1:1%%option to run multiple iterations
                [clusterContainer s accept gammaMatrix] = SplitMergeStep(gammaMatrix, phiVector, s, clusterContainer, k);
                splitMergeCounter                       = splitMergeCounter + [accept, 1];
                %%IF REQUIRED, UPDATE THE GAMMA-RELATED BOOK-KEEPING ARRAYS
                if accept
                    [columnProductMatrix productMatrix sumsOfColumnProducts uniquePhiCoeffVector] = ComputeDerivedGammaArrays(gammaMatrix, phiVector);
                end
            end
        end
    end
    %%----------------------------------------------------------------------
    %% WRITE THE PHIS, MASS PARAMETERS AND CLUSTERIDS TO FILE --------------
    %%----------------------------------------------------------------------
    if(mod(it, samplingFrequency) == 0)
        appendToFile(mcmcFile, phiVector, massParam, s);
        if any(featureSelectionSwitches)
            WriteFeatureParametersToFile(biomarkerFile, clusterContainer);
        end
    end
    %%----------------------------------------------------------------------
    %% IF WE HAVE ANY TIMECOURSE DATA SETS, WE NEED TO DRAW NEW PARAMETERS -
    %%----------------------------------------------------------------------
    if(~isempty(whichTimeCourseDataSets))
        if(mod(it, hyperSamplingFrequency) == 0)
            for k = whichTimeCourseDataSets
                [ nHyperProposals, nHyperAcceptances, clusterContainer] = ...
                    DrawNewHypers(clusterContainer, k, s, dataTypes{k}, fHandles{k}, ...
                    nHyperProposalsMatrix(k,:), nHyperAcceptancesMatrix(k,:));
                nHyperProposalsMatrix(k,:)   = nHyperProposalsMatrix(k,:)   + nHyperProposals;
                nHyperAcceptancesMatrix(k,:) = nHyperAcceptancesMatrix(k,:) + nHyperAcceptances;
            end
        end
    end
    %%----------------------------------------------------------------------
    %% PERIODICALLY SAVE MATLAB VARIABLES TO FILE --------------------------
    %%----------------------------------------------------------------------
    if(mod(it, samplingFrequency) == 0)
        save([outputPath '.mat'], 'vStar', 'massParam', 'phiVector', ...
            'gammaMatrix', 'columnProductMatrix', 'sumsOfColumnProducts', ...
            'uniquePhiCoeffVector', 's', 'clusterContainer', ...
            'productMatrix', 'normalizingConstant')
    end
end

end
%%*****************************************************************************
%%*** END OF MDI.m ************************************************************
%%*****************************************************************************

%%-----------------------------------------------------------------
%%SAMPLE HYPERPARAMETERS FOR TIME COURSE CASE ---------------------
%%-----------------------------------------------------------------
function [ nHyperProposals, nHyperAcceptances, clusterContainer] = ...
    DrawNewHypers(clusterContainer, k, s, dataType, fHandle, nHyperProposals, nHyperAcceptances)

structOfClusters = clusterContainer(k).clusterStruct;
nGenes = clusterContainer(k).nGenes;
clusterIDs = s(:,k);


% Which clusters are currently occupied?
occupiedClusterIDs  = transpose(unique(clusterIDs));


for i = occupiedClusterIDs
    
    for j = 1:3
        nHyperProposals(j) = nHyperProposals(j) + 1;
        currentCluster   = structOfClusters(i);
        currentLogHypers = currentCluster.logHypers;
        currentlogPriorOfLogHypers = currentCluster.logPriorOfLogHypers;
        currentLogMarginalLikelihood = currentCluster.logMarginalLikelihood;
        
        proposedCluster   = currentCluster;
        proposedCluster.covarianceMatrixInverses(1:nGenes) = struct('invertedCovarianceMatrix', [], 'determinant', []);
        proposedLogHypers = currentLogHypers;
        if (j == 1)
            proposedLogHypers(j) = proposedLogHypers(j) + (1.0*randn);
        else
            proposedLogHypers(j) = proposedLogHypers(j) + (0.5*randn);
        end
        proposedCluster.logHypers = proposedLogHypers;
        squaredHypers = currentCluster.squaredHypers;
        squaredHypers(j) = exp(2*proposedLogHypers(j));
        proposedCluster.squaredHypers = squaredHypers;
        l2          = proposedCluster.squaredHypers(1);
        sf2         = proposedCluster.squaredHypers(2);
        timeDiffs   = proposedCluster.timeDiffs;
        if(strcmp(dataType, 'TimeCourseEqualSpacing'))
            proposedCluster.firstRowOfCovarianceMatrix = sf2*exp(-timeDiffs.^2/(2*l2));
        else
            lowerTriangularLogicalMatrix = proposedCluster.lowerTriangularLogicalMatrix;
            covarianceMatrix = sf2*exp(timeDiffs/(2*l2));
            lowerTriangularPart = covarianceMatrix(lowerTriangularLogicalMatrix);
            proposedCluster.lowerTriangularPartOfCovarianceMatrix = lowerTriangularPart;
        end
        
        proposedLogPriorOfLogHypers = currentlogPriorOfLogHypers;
        hyperPriorParameters = proposedCluster.hyperPriorParams;
        proposedLogPriorOfLogHypers(j) = log(normpdf(proposedLogHypers(j), ...
            hyperPriorParameters(j,1), hyperPriorParameters(j,2)  ));
        proposedCluster.logPriorOfLogHypers = proposedLogPriorOfLogHypers;
        proposedCluster = feval(fHandle, proposedCluster, 'invert');
        proposedCluster = feval(fHandle, proposedCluster, 'marginal');
        proposedLogMarginalLikelihood = proposedCluster.logMarginalLikelihood;
        
        logRatio = proposedLogMarginalLikelihood + ...
            sum(proposedLogPriorOfLogHypers) - ...
            currentLogMarginalLikelihood - ...
            sum(currentlogPriorOfLogHypers);
        
        if rand < exp(logRatio)
            structOfClusters(i) = proposedCluster;
            nHyperAcceptances(j) = nHyperAcceptances(j) + 1;
        end
        %disp(structOfClusters(i).logHypers)
    end
    
    
    
    
end
clusterContainer(k).clusterStruct = structOfClusters;
end



%%-----------------------------------------------------------------
%%CALCULATE columnProductIndexMatrix ------------------------------
%%-----------------------------------------------------------------
%%columnProductIndexMatrix -- describes the matrix columns that
%%must be multiplied together in order to find all column prods
function columnProductIndexMatrix = calculateColumnProductIndexMatrix()
global twoPowKminus1 powersOfTwo

columnProductIndexMatrix = zeros(twoPowKminus1,3);
for i = 1:twoPowKminus1   %twoPowKminusKminus1
    remainder = rem(log2(i), 1);
    if remainder == 0
        counter = 0;
        columnProductIndexMatrix(i:end,1) = i;
    else
        counter = counter + 1;
        columnProductIndexMatrix(i,2) = counter;
    end
end
columnProductIndexMatrix(powersOfTwo,:) = [];
columnProductIndexMatrix(:,3) = columnProductIndexMatrix(:,1) + columnProductIndexMatrix(:,2);
%%Example of use:
%%Suppose we want to find the product of the first and third
%%columns of a matrix
%%1) Make a zero vector, v, of length K
%%2) Set v(1) = 1, v(3) = 1
%%3) Reverse v, so that v = v(end:-1:1)
%%4) Convert v from binary to decimal, to get bin2dec(v)
% 5) Look up column bin2dec(v) in columnProductMatrix
end

%%----------------------------------------------------------------------
%%CALCULATE contextInvolvementMatrix -----------------------------------
%%----------------------------------------------------------------------
%%contextInvolvementIndexMatrix -- for each context, tells us
%%whether or not the jth column of columnProductMatrix involves
%%context j
function contextInvolvementIndexMatrix = calculateContextInvolvementIndexMatrix()
global twoPowKminus1 powersOfTwo K
contextInvolvementIndexMatrix = false(K,twoPowKminus1);
for i = 1:K
    currentPowTwo = powersOfTwo(i);
    ending = -1;
    while(ending < twoPowKminus1)
        starting = ending+1+currentPowTwo;
        ending   = starting+currentPowTwo-1;
        contextInvolvementIndexMatrix(i,starting:ending) = true;
    end
end
end


%%-----------------------------------------------------------------
%%CALCULATE phiIndexMatrix -----------------------------------
%%-----------------------------------------------------------------
%%phiIndexMatrix -- describes the order of the phi combinations
%%that is used throughout.  e.g. (K=4) phi_12, phi_13, phi_14,
%%phi_23, phi_24, phi_34.
function phiIndexMatrix = calculatePhiIndexMatrix()
global K
[c r] = meshgrid(1:K);
cUpperTriangle = tril(c,-1);
rUpperTriangle = tril(r,-1);
index1 = cUpperTriangle(~triu(ones(size(cUpperTriangle))))';
index2 = rUpperTriangle(~triu(ones(size(rUpperTriangle))))';
phiIndexMatrix = [index1', index2'];
end


%%-----------------------------------------------------------------
%%CALCULATE coefficientIndexMatrix -----------------------------------
%%-----------------------------------------------------------------
%%coefficientIndexMatrix -- for each phi combination, prescribes which
%%columns of sumOfColumnProducts need to be multiplied together
%%in order to find the normalizing constant
function [ uniqueCoefficientIndexMatrix J] =...
    calculateCoefficientIndexMatrix()
global K twoPowNumberOfPhis twoPowKminus1 powersOfTwo binaryMatrixNumberOfPhis phiIndexMatrix
bigCell = cell(twoPowNumberOfPhis,1);   % We can delete this later, since
% we only need the coefficientIndexMatrix, but for now it helps to see how things are
% working.
coefficientIndexMatrix = false(twoPowNumberOfPhis, twoPowKminus1);
coefficientIndexMatrix(1,powersOfTwo) = true;
for i = 2:twoPowNumberOfPhis
    currentRow = binaryMatrixNumberOfPhis(i,:);
    phiCurrents= phiIndexMatrix(currentRow == 1,:);
    myCell = {};
    counter = 1;
    leftOutContexts = setdiff(1:K, reshape(phiCurrents,1,[]));
    while(~isempty(phiCurrents))
        phiCurrentSet = phiCurrents(1,:);
        phiCurrents(1,:) = [];
        myCell{counter} = phiCurrentSet;
        while(~isempty(phiCurrentSet))
            newPhiSet = [];
            for j = 1:length(phiCurrentSet)
                newInds = union(find(phiCurrents(:,1)==phiCurrentSet(j)),find(phiCurrents(:,2)==phiCurrentSet(j)))';
                if(isempty(newInds))
                else
                    newPhiSet = setdiff([newPhiSet, reshape(phiCurrents(newInds,:),1,[])],phiCurrentSet);
                    phiCurrents(newInds,:) = [];
                end

            end
            
            phiCurrentSet = newPhiSet;
            myCell{counter} = union(myCell{counter}, phiCurrentSet);
        end
        counter = counter + 1;
    end
    bigCell{i} = myCell;
    for j=1:length(myCell)
        v = zeros(1,K);
        v(myCell{j}) = 1;
        v= v(end:-1:1);
        coefficientIndexMatrix(i,bin2dec(num2str(v))) = true;
    end
    coefficientIndexMatrix(i,2.^(leftOutContexts-1)) = true;
end

%%Note that coefficientIndexMatrix has many repeated rows
[uniqueCoefficientIndexMatrix,I,J] = unique(coefficientIndexMatrix, 'rows');
%%note that
%%uniqueCoefficientIndexMatrix(J,:) = coefficientIndexMatrix,
%%coefficientIndexMatrix(I,:)   = uniqueCoefficientIndexMatrix


end


%%-----------------------------------------------------------------
%%CALCULATE contextColumnProductIndexMatrix -----------------------
%%-----------------------------------------------------------------
%%contextColumnProductIndexMatrix -- rows indexed by contexts,
%%columns by unique coefficients.  The ij-entry tells us which
%%column of productColumnMatrix contains the contribution of
%%context i to the jth row of uniqueCoefficientIndexMatrix


function contextColumnProductIndexMatrix = calculateContextColumnProductIndexMatrix()
global uniqueCoefficientIndexMatrix contextInvolvementIndexMatrix K
contextColumnProductIndexMatrix = zeros(K,size(uniqueCoefficientIndexMatrix,1));
for i = 1:K
    currentRow = contextInvolvementIndexMatrix(i,:);
    for j = 1:size(uniqueCoefficientIndexMatrix,1)
        contextColumnProductIndexMatrix(i,j) = find(currentRow.*uniqueCoefficientIndexMatrix(j,:));
    end
end
end





function doNotPertainToContexti = calculateContextPertinentPhis()
global phiIndexMatrix K numberOfPhis
pertainToContexti      = zeros(K, K-1);
doNotPertainToContexti = zeros(K, numberOfPhis - K + 1);
if K == 1
    return
end
for k =1:K
    pertainToContexti(k,:) =...
        union(find(phiIndexMatrix(:,1)==k),find(phiIndexMatrix(:,2)==k));
    setDifference = setdiff(1:numberOfPhis, pertainToContexti(k,:));
    if(~isempty(setDifference))
        doNotPertainToContexti(k,:) = setDifference;
    end
end

end


%%-----------------------------------------------------------------
%%CALCULATE finalIndexMatrix -----------------------
%%-----------------------------------------------------------------
%%finalIndexMatrix -- a simple dictionary to enable us to look
%%up the appropriate phi weighting factors for any
%%configuration of the component labels

function finalIndexMatrix = calculateFinalIndexMatrix()
global K  binaryMatrixKWithoutFirstRow phiIndexMatrix numberOfPhis
finalIndexMatrix = false((2^K)-1, numberOfPhis);

for i = 1:((2^K)-1)
    currentRow = binaryMatrixKWithoutFirstRow(i,:);
    identicalLabels = find(currentRow == 1);
    numIdenticalLabels = length(identicalLabels);
    if(numIdenticalLabels > 1)
        
        l = ones(numIdenticalLabels,1)*identicalLabels;
        r = identicalLabels'*ones(1,numIdenticalLabels);
        L = tril(l,-1);
        R = tril(r,-1);
        phiIndices = [L(L~=0) R(R~=0)];
        finalIndexMatrix(i,ismember(phiIndexMatrix, phiIndices, 'rows')) = true;
        
        
    end
end

end



%% -----------------------------------------------------------------
%%DRAW NEW massParam           ------------------------------
%%-----------------------------------------------------------------
function    massParam = DrawNewMassParameter(massParam, gammaMatrix)
global N K
%%%% Update the mass parameters:
%%We just need a product of the N gammaMatrix distributions
%%associated with each of the N contexts.  We update using a
%%Metropolis-Hastings step:
priorParams = [2,1/4];

%Just iterate through the K mass parameters:
for i = 1:K
    %disp(i)
    currentGammas = gammaMatrix(:,i);
    currentMass   = massParam(i);
    %curLikelihood = prod(gampdf(currentGammas, currentMass/N));
    curLogLikelihood = gamlike([currentMass/N 1], currentGammas);%sum(log(gampdf(currentGammas, currentMass/N)))
    
    %curPrior      = gampdf(currentMass,2,1/4);
    curLogPrior      = gamlike(priorParams, currentMass); %log(gampdf(currentMass,2,1/4));
    
    proposedMass    = currentMass+randn/10;
    if(proposedMass <= 0)
        acceptanceRatio = 0;
    else
        %newLikelihood   = prod(gampdf(currentGammas, proposedMass/N));
        newLogLikelihood = gamlike([proposedMass/N 1], currentGammas);%%sum(log(gampdf(currentGammas, proposedMass/N)));
        
        
        %newPrior        = gampdf(proposedMass,2,1/4);
        newLogPrior        = gamlike(priorParams, proposedMass);%log(gampdf(proposedMass,2,1/4));
        
        %acceptanceRatio = (newLikelihood*newPrior)/(curLikelihood*curPrior);
        acceptanceRatio = exp(-newLogLikelihood-newLogPrior+...
            curLogLikelihood+curLogPrior);
    end
    if(isnan(acceptanceRatio))
        disp(acceptanceRatio)
    end
    
    if(rand < acceptanceRatio)
        massParam(i) = proposedMass;
    end
end

end

%% -----------------------------------------------------------------
%%DRAW NEW vStar PARAMETER           ------------------------------
%%-----------------------------------------------------------------
function vStar = DrawNewVStar(normalizingConstant)
global nGenes
vStar        = gamrnd(nGenes,1/normalizingConstant);
end

%% NOTE:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Recall, if all X_i ~ Exp(b) and all X_i independent, then sum(X_i)
%%  ~ Gam(n, b)  (where i runs from 1 to n)



%% -----------------------------------------------------------------
%%DRAW NEW gammaMatrix PARAMETERS          ------------------------------
%%-----------------------------------------------------------------
function [gammaMatrix columnProductMatrix sumsOfColumnProducts uniquePhiCoeffVector]  = ...
    DrawNewGamma(vStar, massParam, phiVector, gammaMatrix, columnProductMatrix, sumsOfColumnProducts, uniquePhiCoeffVector, s)

global K N J oneMatrix  ...
    contextColumnProductIndexMatrix binaryMatrixNumberOfPhis ...
    powersOfTwo twoPowKminusKminus1 columnProductIndexMatrix

% For each context, we need to keep a tally of the number
% of each of the mixtureIDs
mixtureTallyMatrix = zeros(N,K);
for i = 1:N
    mixtureTallyMatrix(i,:) = histc(s,i);  %Note, we could avoid recalculating this each time, but it justs makes the code more complicated
end

for l = 1:K
    for jl = 1:N
        nl = mixtureTallyMatrix(jl,l);
        
        alteredUniquePhiCoeffVector     = (uniquePhiCoeffVector)./(sumsOfColumnProducts(contextColumnProductIndexMatrix(l,:)));
        
        conditionalUniquePhiCoeffVector = columnProductMatrix(jl,contextColumnProductIndexMatrix(l,:)).*alteredUniquePhiCoeffVector;
        conditionalPhiCoeffVector       = conditionalUniquePhiCoeffVector(J);
        gammaOld = gammaMatrix(jl,l);
        
        meanStar = vStar*sum(conditionalPhiCoeffVector'.*prod(  1+((phiVector(oneMatrix,:)-1).*binaryMatrixNumberOfPhis)  ,2))/gammaOld;
        
        % We then sample the new gamma_{l, jl}
        
        gammaNew = gamrnd(nl + (massParam(l)/N), 1/(meanStar + 1) ) + realmin;
        gammaMatrix(jl,l) = gammaNew;
        
        
        
        % Having sampled this new gammaMatrix, we now need to update the
        % columnProductMatrix and the sumsOfColumnProductMatrix
        
        % %                 % Before the update..
        % %                 oldRowOfColumnProductMatrix = columnProductMatrix(jl,:);
        
        % Now update the appropriate row
        columnProductMatrix(jl,powersOfTwo) = gammaMatrix(jl,:);
        % Calculate column products
        for i = 1:twoPowKminusKminus1
            columnProductMatrix(jl,columnProductIndexMatrix(i,3)) = columnProductMatrix(jl,columnProductIndexMatrix(i,1)).*columnProductMatrix(jl,columnProductIndexMatrix(i,2));
        end
        
        % %                 % After the update...
        % %                 newRowOfColumnProductMatrix = columnProductMatrix(jl,:);
        
        % Update the sumsOfColumnProducts
        % %                 sumsOfColumnProducts = sumsOfColumnProducts -...
        % %                     oldRowOfColumnProductMatrix + newRowOfColumnProductMatrix;
        sumsOfColumnProducts        =  sum(columnProductMatrix);
        
        %sum(columnProductMatrix);
        
        %        disp(sum(abs(sumsOfColumnProducts - sum(columnProductMatrix))))
        
        %         if(K~= 2)
        % Also need to update uniquePhiCoeffVector -- we undo the
        % division that we performed when defining
        % alteredUniquePhiCoeffVector, but use the new column product sums
        uniquePhiCoeffVector = alteredUniquePhiCoeffVector.*sumsOfColumnProducts(contextColumnProductIndexMatrix(l,:));
        %         end
    end
end

% I have commented out the below, as phiCoeffVector and productMatrix are
% recalculated in BlockLabelUpdate anyway, and are not used before then
% % % % Update productMatrix:
% % % phiCoeffVector       = uniquePhiCoeffVector(J);
% % % productMatrix = phiCoeffVector'.*prod(1+...
% % %     ((phiVector(oneMatrix,:)-1).*binaryMatrixNumberOfPhis),2);
end

%% -----------------------------------------------------------------
%%DRAW NEW  phi PARAMETERS -----------------------------------
%%-----------------------------------------------------------------

% Which terms in the normalizing constant involve, say, phi_{1,2} ?  Well,
% it's whichever of the terms has binaryMatrixNumberOfPhis(:,1) == 1.
function [phiVector  productMatrix] = DrawNewPhi(phiVector, vStar, s, productMatrix)
global numberOfPhis phiIndexMatrix binaryMatrixNumberOfPhis


for i = 1:numberOfPhis
    context1 = phiIndexMatrix(i,1);
    context2 = phiIndexMatrix(i,2);
    currentMixtureIDs = [s(:,context1), s(:,context2)];
    phiCurrent   = phiVector(i);
    alphastar    = 1+length(find(currentMixtureIDs(:,1) - currentMixtureIDs(:,2) == 0)); %Number of labels that agree in context1 and context2
    betastar     = vStar*sum(productMatrix(binaryMatrixNumberOfPhis(:,i)))/phiCurrent;
    phiNew       = gamrnd(alphastar, 1 / (betastar+0.2));   %%We assume a Ga(1, 0.2) prior for phi
    phiVector(i) = phiNew;
    % We need to update productMatrix
    productMatrix(binaryMatrixNumberOfPhis(:,i),:) = productMatrix(binaryMatrixNumberOfPhis(:,i),:)*phiNew/phiCurrent;
end
end



%%-----------------------------------------------------------------
%% BLOCK UPDATE COMPONENT LABELS -----------------------------------
%%-----------------------------------------------------------------
%%We allow component labels to be swapped (within contexts).  For
%%example, we might swap around labels 1 and 2 in context 1.
function [s clusterContainer gammaMatrix columnProductMatrix productMatrix sumsOfColumnProducts uniquePhiCoeffVector] = ...
    BlockLabelUpdate(s, clusterContainer, phiVector, gammaMatrix)

global K N finalIndexMatrix twoPowKminus1 doNotPertainToContexti nUniqueCoefficients oneMatrix ...
    uniqueCoefficientIndexMatrix J binaryMatrixNumberOfPhis powersOfTwo columnProductIndexMatrix twoPowKminusKminus1
% Now allow us to propose a change which will swap the labels of two
% components in context j
% e.g. if we have
% Context 1: 1 1 1 2 2 2 3 3
% Context 2: 2 2 2 1 1 1 3 3
% then this might swap the 1s and 2s in Context 1
for j = 1:K
    clustersForThisContext    = clusterContainer(j).clusterStruct;
    notPertinentInThisContext = doNotPertainToContexti(j,:);
    for i = 1:N
        % Randomly choose any component that is not component i
        newpos = ceil(rand * (N - 1));
        if ( newpos >= i )
            newpos = newpos + 1;
        end
        %We propose to swap the labels of components i and newpos
        originalInds       = s(:,j) == i;
        originalIndsNewpos = s(:,j) == newpos;
        %Calculate likelihood before the swap:
        labels1 = s(originalInds,:); % This returns all rows of s where the jth column contains i.
        labels2 = s(originalIndsNewpos,:); % This returns all rows of s where the jth column contains newpos.
        %Determine the phi dependence coefficients that we need:
        agreeingLabels      = [(labels1 == i); (labels2 == newpos)];
        binInd              = agreeingLabels*(2.^(size(agreeingLabels,2)-1:-1:0))';  % = bin2dec(num2str(agreeingLabels))
        phiInds             = finalIndexMatrix(binInd,:);
        phiMat              = ones(size(agreeingLabels,1),1)*phiVector;
        logPhiSum           = sum(log(1+phiMat(phiInds)));  %Note that the "1+" converts from phis to rhos
        swappedLabels1      = labels1;
        swappedLabels2      = labels2;
        swappedLabels1(:,j) = newpos;
        swappedLabels2(:,j) = i;
        %Determine the phi dependence coefficients that we need:
        agreeingLabelsSwap = [(swappedLabels1 == newpos); (swappedLabels2 == i)];
        binInd             = agreeingLabelsSwap*(2.^(size(agreeingLabelsSwap,2)-1:-1:0))';% = bin2dec(num2str(agreeingLabelsSwap));
        phiIndsSwap        = finalIndexMatrix(binInd,:);
        phiIndsSwap(:,notPertinentInThisContext) = false;
        phiMatSwap         = ones(size(agreeingLabelsSwap,1),1)*phiVector;
        logPhiSumSwap      = sum(log(1+phiMatSwap(phiIndsSwap)));  %Note that the "1+" converts from phis to rhos
        logaccept          = logPhiSumSwap - logPhiSum;
        
        accept = 1;
        if ( logaccept < 0 )
            accept = exp(logaccept);
        end
        
        if ( rand < accept )
            % swap the labels:
            savedCluster                   = clustersForThisContext(i);
            clustersForThisContext(i)      = clustersForThisContext(newpos);
            clustersForThisContext(newpos) = savedCluster;
            s(originalInds,j)              = newpos;
            s(originalIndsNewpos,j)        = i;
            %We also need to swap the appropriate elements of gammaMatrix:
            savedGamma            = gammaMatrix(i,j);
            gammaMatrix(i,j)      = gammaMatrix(newpos, j);
            gammaMatrix(newpos,j) = savedGamma;
            %This has an impact on columnProductMatrix too, so recalculate.
            %Might as well do this once, at the end.
            %%
            %This, in turn, has an impact on sumsOfColumnProducts, which
            %then has an impact upon productMatrix... but we'll just
            %recalculate these at the end too, rather than doing clever
            %indexing
            %%
            %Note that the commented code below is a faster way to update,
            %but has a tendency to create numerical errors (leading to
            %NaNs)
            % %             columnProductMatrix(i, contextInvolvementIndexMatrix(j,:)) ...
            % %                 = columnProductMatrix(i, contextInvolvementIndexMatrix(j,:))/savedGamma;
            % %             columnProductMatrix(i, contextInvolvementIndexMatrix(j,:)) ...
            % %                 = columnProductMatrix(i, contextInvolvementIndexMatrix(j,:))*gammaMatrix(i,j);
            % %
            % %             columnProductMatrix(newpos, contextInvolvementIndexMatrix(j,:)) ...
            % %                 = columnProductMatrix(newpos, contextInvolvementIndexMatrix(j,:))/gammaMatrix(i,j);
            % %             columnProductMatrix(newpos, contextInvolvementIndexMatrix(j,:)) ...
            % %                 = columnProductMatrix(newpos, contextInvolvementIndexMatrix(j,:))*gammaMatrix(newpos,j);
            
        end
    end
    clusterContainer(j).clusterStruct = clustersForThisContext;
end
%%RECOMPUTE VARIOUS BOOK-KEEPING ARRAYS
[columnProductMatrix productMatrix sumsOfColumnProducts uniquePhiCoeffVector] = ComputeDerivedGammaArrays(gammaMatrix, phiVector);

end
%%-----------------------------------------------------------------
%% DRAW NEW ITEM LABELS --------------------------------------------
%%-----------------------------------------------------------------
function [clusterContainer s] = DrawNewItemLabel(gammaMatrix, phiVector, s, clusterContainer)
%%GLOBAL VARIABLES REQUIRED BY THIS FUNCTION
global K nGenes N phiIndexMatrix finalIndexMatrix fHandles ...
    doNotPertainToContexti allLabelsMatrix timeCourseSwitches ...
    gaussianSwitches poissonSwitches nbSwitches massParam
allComponents = 1:N;
% Start with a fixed number of particles allow specification later
numbofparts = 1;
logweight = zeros(1, numbofparts);
partstar = zeros(1, numbofparts);

%%LOOP OVER ALL 'K' DATA TYPES
for j = 1:K
    if K > 1
        phiIndicesForConditional = union(find(phiIndexMatrix(:,1)== j), find(phiIndexMatrix(:,2)== j));
    else
        phiIndicesForConditional = [];
    end
    notPertinentInThisContext = doNotPertainToContexti(j,:);
    clustersForThisContext    = clusterContainer(j).clusterStruct;
    dataForThisContext        = clusterContainer(j).data;
    fHandle = fHandles{j};
    
    % Declare some variables for the particle filter
    M =massParam(1, j);
    a = 0.5;
    mu = 0;
    sigmasq = 1;
    % Use nGenes for n
    % Make these greater dimension to allow for the multiple datasets; not sure yet
    sumy = cell(1, numbofparts);
    sumysq = cell(1, numbofparts);
    nj = cell(1, numbofparts);
    sumv = zeros(numbofparts, 1);
    % Need to declare s...?
    sstar = zeros(numbofparts, nGenes);
    %%LOOP OVER ALL GENES IN EACH TYPE
    for i = 1:nGenes
        % disp(['i = ' num2str(i)]);
        % This is where the particle filter needs to come in
        % Loop over each of the particles
        % First step is draw the s values
        % The allocation variables conditional on y[1:t] and s*(particle)[1:(t-1)]
        % q(k) poportional to m(i){k, t-1}k{*}{k}(yt|s*(i){1:(t-1)) if k between 1 and Ki{t-1}
        % else if k = Ki{t - 1} + 1 M k*new(yt)
        % Need to update fprob using the particle filter
        %%FIND USEFUL VALUES
        oldLabel = s(i,j);
        % disp(['oldlabel' num2str(oldLabel)]);
        prob     = gammaMatrix(:, j);%%conditional prior
        %%We need to find the labels of gene i in the other contexts:
        %%We consider all possibilities for the label of gene i
        %%in context j:
        labelsAcrossAllContexts            = s(i,:);
        labelsAcrossAllContextsMatrix      = ones(N,1)*labelsAcrossAllContexts;
        labelsAcrossAllContextsMatrix(:,j) = 1:N;
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
        %%NOW FIND THE "LIKELIHOOD" TERM
        proposedClustersForThisContext = clustersForThisContext;
        occupiedClusters               = transpose(unique(s(:,j)));
        %%ADD/REMOVE ITEMS FOR PROPOSED MOVES
        %%the removal means our proposed marginal likelihood for the
        %current cluster will be the inverse of what we want, but this
        %is corrected for below
        for ind = occupiedClusters
        if(ind == oldLabel), proposedClustersForThisContext = AddRemoveItem('removeGene', proposedClustersForThisContext, ind, i, dataForThisContext, j);
        else                 proposedClustersForThisContext = AddRemoveItem('addGene',    proposedClustersForThisContext, ind, i, dataForThisContext, j);, end
        % We have inverted a covariance matrix, so should not let this go to waste...
          if(timeCourseSwitches(j))
              nGenesInProposedCluster = proposedClustersForThisContext(ind).nGenes;
              if(nGenesInProposedCluster > 0)
                  clustersForThisContext(ind).covarianceMatrixInverses(nGenesInProposedCluster) = ...
                      proposedClustersForThisContext(ind).covarianceMatrixInverses(nGenesInProposedCluster);
              end
          end
        end
        %%WE ALSO NEED THE POSSIBILITY OF ADDING THE ITEM TO AN UNOCCUPIED CLUSTER
        unoccupiedClusters     = setdiff(allComponents,occupiedClusters);
        firstUnoccupiedCluster = unoccupiedClusters(1);
        if(timeCourseSwitches(j))
            proposedClustersForThisContext(firstUnoccupiedCluster).covarianceMatrixInverses(1:nGenes) = struct('invertedCovarianceMatrix', [], 'determinant', []);
        end
        proposedClustersForThisContext                     = AddRemoveItem('addGene', proposedClustersForThisContext, firstUnoccupiedCluster, i, dataForThisContext, j);
        proposedClustersForThisContext(unoccupiedClusters) = proposedClustersForThisContext(firstUnoccupiedCluster);
        % Loop over the particles
        for it = 1:numbofparts
            if(i == 1)
                sumy{1, it} = dataForThisContext(i);
                sumysq{1, it} = dataForThisContext(i) ^ 2;
                nj{1, it} = 1;
                logweight(it) = 0;
                sstar(it, i) = i;
            else
                prob = [nj{1, it} M] ./ (sum(nj{1, it}) + M);
                % disp(['m' num2str(M)]);
                % disp(['prob = ' num2str(prob)]);
                % Might want to draw these mu star values from the previously allocated clustering?
                mustar = [(mu / (1 - a) + sumy{1, it} / a) ./ (nj{1, it} / a + 1 / (1 - a)) mu];
                varstar = [sigmasq ./ (nj{1, it} / a + 1 / (1 - a)) sigmasq*(1-a)];
                logprob = - 0.5 * (dataForThisContext(i) - mustar).^2 ./ (a * sigmasq + varstar) - 0.5 * log(a * sigmasq + varstar);
                % disp(['logprob' num2str(logprob)]);
                fprob = cumsum(prob .* exp(logprob - max(logprob)));
                logweight(it) = log(fprob(end)) + max(logprob);
                fprob = fprob / fprob(end);
                u1 = rand;
                sstar(it, i) = 1;
                while ( fprob(sstar(it, i)) < u1 )
                    sstar(it, i) = sstar(it, i) + 1;
                end
                if ( sstar(it, i) == length(nj{1, it}) + 1 )
                    nj{1, it} = [nj{1, it} 1];
                    sumy{1, it} = [sumy{1, it} dataForThisContext(i)];
                    sumysq{1, it} = [sumysq{1, it} dataForThisContext(i)^2];
                else
                    nj{1, it}(sstar(it, i)) = nj{1, it}(sstar(it, i)) + 1;
                    sumy{1, it}(sstar(it, i)) = sumy{1, it}(sstar(it, i)) + dataForThisContext(i);
                    sumysq{1, it}(sstar(it, i)) = sumysq{1, it}(sstar(it, i)) + dataForThisContext(i)^2;
                end
            end
        end
        fprob = cumsum(exp(logweight - max(logweight)));
        fprob = fprob / fprob(end);
        % disp(['fprob = ' num2str(fprob)]);
        u1 = rand / numbofparts;
        m = 1;
        it = 1;
        while (m <= numbofparts )
            while ( (u1 < fprob(m)) && (m <= numbofparts) )
                partstar(it) = m;
                u1 = u1 + 1 / numbofparts;
                it = it + 1;
            end
            m = m + 1;
        end
        sumy = sumy(1, partstar);
        sumysq = sumysq(1, partstar);
        nj = nj(1, partstar);
        sumv = sumv(partstar);
        % Here we need to pick out a particular particle surely?
        % But then future iterations of the particle filter will be based only on a single particle
        % Then perhaps the algorithm should only update the clusterstruct after the particles have been fully drawn?   
        s(1:i, j) = transpose(sstar(partstar, 1:i));


%{
    %%COMPUTE THE MARGINAL LIKELIHOOD RATIOS
    oldLogMarginalLikelihoods         = [clustersForThisContext.logMarginalLikelihood];
    newLogMarginalLikelihoods         = [proposedClustersForThisContext.logMarginalLikelihood];
    marginalLikelihoodRatio           = newLogMarginalLikelihoods - oldLogMarginalLikelihoods;
    marginalLikelihoodRatio(oldLabel) = -marginalLikelihoodRatio(oldLabel);%%correction for the add/remove
    %%REMOVE A CONSTANT FACTOR BEFORE EXPONENTIATING
    %%normalising here prevents numerical instability
    marginalLikelihoodRatio = marginalLikelihoodRatio - max(marginalLikelihoodRatio);
    prob                    = prob.*exp(marginalLikelihoodRatio');
    fprob        	           = cumsum(prob);
    fprob                   = fprob / fprob(end);  %Find the (normalised) cumulative distribution
    
    u       = rand;
    s(i, j) = 1;
    while ( u > fprob(s(i, j)) )
        s(i, j) = s(i, j) + 1;   % Allocate the gene to a new component, by sampling according to the probabilities just calculated
    end
    % We now need to update the
    % the component from which item i was removed, and for
    % the one to which it was added (unless we have put the
    % item back in the same cluster):
%}
    %%WE ALSO NEED THE POSSIBILITY OF ADDING THE ITEM TO AN UNOCCUPIED CLUSTER
    unoccupiedClusters     = setdiff(allComponents,occupiedClusters);
    firstUnoccupiedCluster = unoccupiedClusters(1);
      
    ind = s(i, j);
    if(s(i,j) ~= oldLabel)
        clustersForThisContext(oldLabel) = proposedClustersForThisContext(oldLabel);
        clustersForThisContext(s(i,j))   = proposedClustersForThisContext(s(i,j));
    end
end
clusterContainer(j).clusterStruct = clustersForThisContext;
end

end
%%----------------------------------------------------------------------
%% INITIALISE THE MCMC FILE --------------------------------------------
%%----------------------------------------------------------------------
%%> 2) a CSV file that contains the sampled values. Let K=nDatasets. The first K
%%> columns provide the sampled values for the K mass parameters associated with
%%> each Dirichlet prior. The next K(K-1)/2 columns provide the "phi" (dataset
%%> association) parameters. These are in the order phi_{1,2}, phi_{1,3}, ...,
%%> phi_{1,K}, phi_{2,3}, phi_{2, 4}, ... and so on. The next p columns provide
%%> the component allocation labels for the p genes in Dataset 1. The next p
%%> columns contain the component allocation labels for the p genes in Dataset 2.
%> ... and so on.
function writeNewFile(outputFile, phiVector, massParam, s, itemNames)
%%FIND USEFUL VALUES
nDataTypes   = length(massParam);
nDataItems   = length(itemNames);
outputString = [];
%%HEADER:  MASS PARAMETERS
for i=1:nDataTypes
    outputString = [outputString, 'MassParameter_', num2str(i), ','];
end
%%HEADER:  PHI PARAMETERS
for i=1:(nDataTypes-1)
    for j=(i+1):nDataTypes
        outputString = [outputString, 'Phi_', num2str(i), num2str(j), ','];
    end
end
%%HEADER: COMPONENT ALLOCATION LABELS
for i=1:nDataTypes
    for j=1:nDataItems
        outputString = [outputString, 'Dataset', num2str(i), '_', itemNames{j}, ','];
    end
end
outputString = outputString(1:(end-1));%%clip the final comma
outputString = [outputString '\n'];
%%WRITE HEADER TO FILE
fid = fopen(outputFile, 'wt');
fprintf(fid, outputString);
fclose(fid);
%%WRITE THE INITIAL STATE TO FILE
appendToFile(outputFile, phiVector, massParam, s);
end
%%----------------------------------------------------------------------
%% APPEND PARAMETERS TO THE MCMC FILE ----------------------------------
%%----------------------------------------------------------------------
function appendToFile(outputFile, phiVector, massParam, s)
outputLine = [massParam, phiVector, s(:)'];
dlmwrite(outputFile, outputLine, '-append');
end
%%----------------------------------------------------------------------
%% INITIALISE THE OUTPUT FILE CONTAINING FEATURE PARAMETERS ------------
%%----------------------------------------------------------------------
function InitialiseBiomarkerFile(outputFile, clusterContainer)
%%CONSTRUCT OUTPUT HEADER
nDataTypes   = length(clusterContainer);
outputString = [];
for i=1:nDataTypes
    dataName     = ['DataType', num2str(i), '_'];
    featureNames = clusterContainer(i).featureNames;
    featureNames = strcat(dataName, featureNames, ', ');
    newString    = [featureNames{:}];
    outputString = [outputString, newString];
end
outputString = outputString(1:(end-1));%%clip the final comma
outputString = [outputString '\n'];
%%WRITE HEADER TO FILE
fid = fopen(outputFile, 'wt');
fprintf(fid, outputString);
fclose(fid);
%%WRITE THE INITIAL STATE TO FILE
WriteFeatureParametersToFile(outputFile, clusterContainer);
end
%%----------------------------------------------------------------------
%% WRITE THE FEATURE PARAMETERS OUT TO FILE ----------------------------
%%----------------------------------------------------------------------
function WriteFeatureParametersToFile(outputFile, clusterContainer)
nDataTypes = length(clusterContainer);
outputLine = [];
for i=1:nDataTypes
    outputLine = [outputLine, clusterContainer(i).clusterStruct(1).featureParameters];
end
dlmwrite(outputFile, outputLine, '-append')
end
%%----------------------------------------------------------------------
%% ADD/REMOVE ITEM FROM A GIVEN CLUSTER --------------------------------
%%----------------------------------------------------------------------
function clusterStruct = AddRemoveItem(whichOperation, clusterStruct, whichID, index, data, whichType)
%%GLOBAL VARIABLES WE MAY NEED
global fHandles gaussianSwitches poissonSwitches nbSwitches
%%EXTRACT THE REQUIRED DATA
if gaussianSwitches(whichType) | poissonSwitches(whichType) | nbSwitches(whichType)
    % PDWK 20120120
    % This is slightly hacky, since - in this case -
    % "dataForCurrentGene" is not just the data for the
    % current gene, but actually the whole dataset for the
    % context!  However, this is the easiest way to make
    % the changes required for the Gaussian data type
    currentData = data;
else
    currentData = data(index,:);
end
%%ADD OR REMOVE THE ITEM, AS REQUIRED
clusterStruct(whichID) = feval(fHandles{whichType}, clusterStruct(whichID), whichOperation, index, currentData);
end
%%----------------------------------------------------------------------
%% MAKE A SPLIT-MERGE METROPOLIS-HASTINGS STEP -------------------------
%%----------------------------------------------------------------------
%%this is a Metropolis-Hastings MCMC step, based on the
%%split-merge algorithm of Dahl (2005):
%%'Sequentially-allocated merge-split sampler for conjugate and
%nonconjugate Dirichlet process mixture models'
%%
%%We expect these steps to have a low acceptance rate.  They are
%included for the occasional cases where the broad structure of the
%clustering partition is in a sub-optimal region of parameter
%space.  In this case, an occasional "kick" from a split-merge step
%can be very helpful.
%%
%%NOTE: for now, only do this for the conjugate data models
%%the non-conjugate case is a bit more complicated
%%
%%NOTE: We may have updated the gammaMatrix in this function.  If
%so, we will need to update the various associated book-keeping
%arrays.  For convenience, this is NOT done in this function.
%%
function [clusterContainer s accept gammaMatrix] = SplitMergeStep(gammaMatrix, phiVector, s, clusterContainer, whichType);
%%GLOBAL VARIABLES WE MAY NEED
global phiIndexMatrix fHandles
%%FIND USEFUL VALUES
clusterIDs      = s(:, whichType);
nDataItems      = size(s, 1);
nDataTypes      = size(s, 2);
data            = clusterContainer(whichType).data;
currentStruct   = clusterContainer(whichType).clusterStruct;
proposedStruct  = clusterContainer(whichType).clusterStruct;
%%RANDOMLY SELECT TWO ITEMS TO SEED THE SPLIT-MERGE STEP
index_1     = randi(nDataItems);
working     = 1:nDataItems;
working     = working(working~=index_1);
index_2     = working(randi(nDataItems-1));
clusterID_1 = s(index_1, whichType);
clusterID_2 = s(index_2, whichType);
splitSwitch = clusterID_1==clusterID_2;
%%GENERATE A PROPOSAL STEP
if splitSwitch
    [index, newIDs, proposalLogProb, proposedGammas, moveID] = SplitProposal(proposedStruct, clusterIDs, index_1, index_2, data, whichType, gammaMatrix);
else
    [index, newIDs, proposalLogProb, proposedGammas] = MergeProposal(proposedStruct, clusterIDs, index_1, index_2, data, whichType, gammaMatrix);
end
oldIDs = s(index, whichType);
%%MOVE THE REQUIRED ITEMS TO THEIR NEW CLUSTER
%%obv only need to do this is the ID has changed!
for i=1:length(index)
    if newIDs(i)~=oldIDs(i)
        proposedStruct = AddRemoveItem('removeGene', proposedStruct, oldIDs(i), index(i), data, whichType);
        proposedStruct = AddRemoveItem('addGene',    proposedStruct, newIDs(i), index(i), data, whichType);
    end
end
%%THE LOG MARIGNAL LIKELIHOODS
%%we only need the clusters that are changed in this step
if splitSwitch
    logPosteriorRatio = proposedStruct(clusterID_1).logMarginalLikelihood + ...
        proposedStruct(moveID).logMarginalLikelihood      - ...
        currentStruct(clusterID_1).logMarginalLikelihood  - ...
        currentStruct(moveID).logMarginalLikelihood;
else
    logPosteriorRatio = proposedStruct(clusterID_1).logMarginalLikelihood + ...
        proposedStruct(clusterID_2).logMarginalLikelihood - ...
        currentStruct(clusterID_1).logMarginalLikelihood  - ...
        currentStruct(clusterID_2).logMarginalLikelihood;
end
%%EFFECT OF THE MIXTURE WEIGHTS (GAMMA)
%%assume that the sum of gammas hasn't changed...
gamma_old         = gammaMatrix(oldIDs, whichType);
gamma_new         = proposedGammas(newIDs, whichType);
logPosteriorRatio = logPosteriorRatio + sum(log(gamma_new ./ gamma_old));
%%EFFECT OF THE FUSION WEIGHTS (PHI)
%%the moved items will may match clusterIDs in other data types
phiIndex                = 0;
s_new                   = s;
s_new(index, whichType) = newIDs;
for i=1:(nDataTypes-1)
    for j=(i+1):nDataTypes
        phiIndex          = phiIndex+1;
        nPhiCounts        = sum(s_new(:,i)==s_new(:,j));
        nPhiCounts        = nPhiCounts - sum(s(:,i)==s(:,j));
        logPosteriorRatio = logPosteriorRatio + nPhiCounts*log(1+phiVector(phiIndex));
    end
end
%%DECIDE WHETHER TO MAKE THE M-H STEP
logPosteriorRatio = logPosteriorRatio + proposalLogProb;
acceptanceProb    = exp(logPosteriorRatio);
accept            = acceptanceProb > unifrnd(0, 1);
%%IF REQUIRED, UPDATE THE CLUSTER CONTAINERS
if accept
    s                                         = s;
    s(index, whichType)                       = newIDs;
    clusterContainer(whichType).clusterStruct = proposedStruct;
    gammaMatrix                               = proposedGammas;
end

end
%%----------------------------------------------------------------------
%% SPLIT PROPOSAL STEP -------------------------------------------------
%%----------------------------------------------------------------------
%%in this case, both item1 and item2 are initially in the same cluster
%%The proposed step is to randomly assign all the items in this
%cluster to two new clusters, defined by the initial seed items
%%The first of these new clusters will retain the old ID label; the
%new one will be drawn at random
function [index, newIDs, proposalLogProb, proposedGammas, proposalID_2] = SplitProposal(clusterStruct, clusterIDs, index_1, index_2, data, whichType, gammaMatrix);
%%FIND USEFUL VALUES
nMixtures       = length(clusterStruct);
initialID       = clusterIDs(index_1);
proposalID_1    = initialID;
proposalLogProb = 0;
%%DRAW A NEW CLUSTER TO WHICH TO POPULATE
proposalID_2 = find(ismember(1:nMixtures, unique(clusterIDs))==0);
proposalID_2 = proposalID_2(1);
%%FIND THE RELEVANT ITEMS; RANDOMISE THEIR ORDER
%%exclude the seed items from these arrays initially, becasue we're
%automatically starting with these.  Add them after randomly
%shuffling the order
index            = find(clusterIDs==initialID);
index            = index(index~=index_1);
index            = index(index~=index_2);
newIDs           = clusterIDs(index);
nMoveItems       = length(index);
[dum, sortIndex] = sort(rand(1, nMoveItems));
index            = index(sortIndex);
newIDs           = newIDs(sortIndex);
index            = [index_1,   index_2,   index'];
newIDs           = [initialID, initialID, newIDs'];
nMoveItems       = length(index);
%%REMOVE THE CURRENT ITEMS FROM THE CLUSTER STRUCTURE
%%the first item stays where it is, by construction
%%the second item must move to the new cluster
for i=2:nMoveItems
    clusterStruct = AddRemoveItem('removeGene', clusterStruct, newIDs(i), index(i), data, whichType);
end
clusterStruct = AddRemoveItem('addGene', clusterStruct, proposalID_2, index(2), data, whichType);
newIDs(2)     = proposalID_2;
%%PROBABILISTICALLY ASSIGN EACH OTHER ITEM TO ONE OF THE TWO SPLIT CLUSTERS
for i=3:nMoveItems
    %%COMPUTE THE PROBABILITY OF ASSIGNING THE ITEM TO EACH CLUSTER
    dum1 = AddRemoveItem('addGene', clusterStruct, proposalID_1, index(i), data, whichType);
    dum2 = AddRemoveItem('addGene', clusterStruct, proposalID_2, index(i), data, whichType);
    %%COMPUTE THE NORMALISED MOVE PROBABILITY
    stayLogProb = dum1(proposalID_1).logMarginalLikelihood - clusterStruct(proposalID_1).logMarginalLikelihood;
    moveLogProb = dum2(proposalID_2).logMarginalLikelihood - clusterStruct(proposalID_2).logMarginalLikelihood;
    normalise   = max(stayLogProb, moveLogProb);
    stayLogProb = stayLogProb - normalise;
    moveLogProb = moveLogProb - normalise;
    normalise2  = log(exp(stayLogProb) + exp(moveLogProb));
    stayLogProb = stayLogProb - normalise2;
    moveLogProb = moveLogProb - normalise2;
    moveProb    = exp(moveLogProb);
    %%RANDOMLY ASSIGN THE ITEM TO ONE OF THE CLUSTERS; UPDATE PROPOSAL PROB
    if rand()<moveProb
        clusterStruct   = dum2;
        newIDs(i)       = proposalID_2;
        proposalLogProb = proposalLogProb - moveLogProb;
    else
        clusterStruct   = dum1;
        proposalLogProb = proposalLogProb - stayLogProb;
    end
end
%%ALSO PROPOSE A STEP IN THE GAMMA WEIGHTS
%%the proposalLogProb will also need updating
%%as will some book-keeping arrays, I think...
%%This step simple redistributes the gamm probability masses of the
%two clusters according to a beta distribution, using the nItems in
%each cluster
proposedGammas                          = gammaMatrix;
gamma_1                                 = proposedGammas(proposalID_1, whichType);
gamma_2                                 = proposedGammas(proposalID_2, whichType);
totalGamma                              = gamma_1 + gamma_2;
%%BETA DISTRIBUTED PROPOSAL
nItems_1                                = clusterStruct(proposalID_1).nGenes;
nItems_2                                = clusterStruct(proposalID_2).nGenes;
nTotalItems                             = nItems_1 + nItems_2;
randomSplit                             = betarnd(1+nItems_1, 1+nItems_2);
proposedGammas(proposalID_1, whichType) = randomSplit*totalGamma;
proposedGammas(proposalID_2, whichType) = (1-randomSplit)*totalGamma;
%%also need to do the proposal distribution...
newLogProb      = log(betapdf(gamma_1/totalGamma, 1+nTotalItems, 1));
newLogProb      = newLogProb - log(betapdf(randomSplit,   1+nItems_1, 1+nItems_2));
proposalLogProb = proposalLogProb + newLogProb;
end
%%----------------------------------------------------------------------
%% MERGE PROPOSAL STEP -------------------------------------------------
%%----------------------------------------------------------------------
%%in this case, item1 and item2 are initially in different clusters
%%The proposed step is to merge all items into a single cluster,
%using the label of item1
function [index, newIDs, proposalLogProb, proposedGammas] = MergeProposal(clusterStruct, clusterIDs, index_1, index_2, data, whichType, gammaMatrix);
%%FIND USEFUL VALUES
nMixtures       = length(clusterStruct);
initialID_1     = clusterIDs(index_1);
initialID_2     = clusterIDs(index_2);
proposalID      = initialID_1;
proposalLogProb = 0;
%%FIND THE RELEVANT ITEMS; RANDOMISE THEIR ORDER
%%exclude the seed items from these arrays initially, becasue we're
%automatically starting with these.  Add them after randomly
%shuffling the order
index            = find(clusterIDs==initialID_1 | clusterIDs==initialID_2);
index            = index(index~=index_1);
index            = index(index~=index_2);
newIDs           = clusterIDs(index);
nMoveItems       = length(index);
[dum, sortIndex] = sort(rand(1, nMoveItems));
index            = index(sortIndex);
newIDs           = newIDs(sortIndex);
index            = [index_1,     index_2,     index'];
initialIDs       = [initialID_1, initialID_2, newIDs'];
newIDs           = initialIDs;
nMoveItems       = length(index);
%%REMOVE THE CURRENT ITEMS FROM THE CLUSTER STRUCTURE
%%the first item stays where it is, by construction
%%the second item must move to the new cluster
for i=2:nMoveItems
    clusterStruct = AddRemoveItem('removeGene', clusterStruct, newIDs(i), index(i), data, whichType);
end
clusterStruct = AddRemoveItem('addGene', clusterStruct, proposalID, index(2), data, whichType);
newIDs(2)     = proposalID;
%%TO ASSESS THE PROPOSAL PROBABILITY, WE NEED TO CONSIDER THE REVERSE STEP
for i=3:nMoveItems
    %%COMPUTE THE PROBABILITY OF ASSIGNING THE ITEM TO EACH CLUSTER
    dum1 = AddRemoveItem('addGene', clusterStruct, initialID_1, index(i), data, whichType);
    dum2 = AddRemoveItem('addGene', clusterStruct, initialID_2, index(i), data, whichType);
    %%COMPUTE THE NORMALISED MOVE PROBABILITY
    stayLogProb = dum1(initialID_1).logMarginalLikelihood - clusterStruct(initialID_1).logMarginalLikelihood;
    moveLogProb = dum2(initialID_2).logMarginalLikelihood - clusterStruct(initialID_2).logMarginalLikelihood;
    normalise   = max(stayLogProb, moveLogProb);
    stayLogProb = stayLogProb - normalise;
    moveLogProb = moveLogProb - normalise;
    normalise2  = log(exp(stayLogProb) + exp(moveLogProb));
    stayLogProb = stayLogProb - normalise2;
    moveLogProb = moveLogProb - normalise2;
    moveProb    = exp(moveLogProb);
    %%UPDATE DEPENDING ON WHERE THE ITEM COMES FROM
    if initialIDs(i)~=proposalID
        proposalLogProb = proposalLogProb + moveLogProb;
        clusterStruct   = dum2;
    else
        proposalLogProb = proposalLogProb + stayLogProb;
        clusterStruct   = dum1;
    end
    %%UPDATE THE ACTUAL MOVE
    newIDs(i)     = proposalID;
end
%%ALSO PROPOSE A STEP IN THE GAMMA WEIGHTS
%%the proposalLogProb will also need updating
%%as will some book-keeping arrays, I think...
%%This step simple redistributes the gamm probability masses of the
%two clusters according to a beta distribution, using the nItems in
%each cluster
proposedGammas                         = gammaMatrix;
gamma_1                                = proposedGammas(initialID_1, whichType);
gamma_2                                = proposedGammas(initialID_2, whichType);
totalGamma                             = gamma_1 + gamma_2;
%%BETA DISTRIBUTED PROPOSAL
nItems_1                               = clusterStruct(initialID_1).nGenes;
nItems_2                               = clusterStruct(initialID_2).nGenes;
nTotalItems                            = nItems_1 + nItems_2;
randomSplit                            = betarnd(1+nTotalItems, 1);
proposedGammas(initialID_1, whichType) = randomSplit*totalGamma;
proposedGammas(initialID_2, whichType) = (1-randomSplit)*totalGamma;
%%also need to do the proposal distribution...
newLogProb      = log(betapdf(gamma_1/totalGamma,       1+nItems_1, 1+nItems_2));
newLogProb      = newLogProb - log(betapdf(randomSplit, 1+nTotalItems, 1));
proposalLogProb = proposalLogProb + newLogProb;
end
%%----------------------------------------------------------------------
%% RECOMPUTE THE DERIVED GAMMA ARRAYS ----------------------------------
%%----------------------------------------------------------------------
%%This function computes a set of book-keeping arrays derived from
%gammaMatrix.  These must be recomputed whenever the gammas are resampled.
function [columnProductMatrix productMatrix sumsOfColumnProducts uniquePhiCoeffVector] = ComputeDerivedGammaArrays(gammaMatrix, phiVector)
%%GLOBAL VARIABLES WE REQUIRE
global K N finalIndexMatrix twoPowKminus1 doNotPertainToContexti nUniqueCoefficients oneMatrix ...
    uniqueCoefficientIndexMatrix J binaryMatrixNumberOfPhis powersOfTwo columnProductIndexMatrix twoPowKminusKminus1
%%FIND USEFUL VALUES
columnProductMatrix                = zeros(N,twoPowKminus1);
columnProductMatrix(:,powersOfTwo) = gammaMatrix;
%%COMPUTE THE COLUMN PRODUCTS
for i = 1:twoPowKminusKminus1
    columnProductMatrix(:,columnProductIndexMatrix(i,3)) = columnProductMatrix(:,columnProductIndexMatrix(i,1)).*columnProductMatrix(:,columnProductIndexMatrix(i,2));
end
sumsOfColumnProducts = sum(columnProductMatrix);
%%ERROR-CHECKING
%usually because the gammas have become too small, and caused numerical issues!
if(sum(isnan(columnProductMatrix)))
    error('debug')
end
%%ALSO COMPUTE OTHER ARRAYS OF INTEREST
uniquePhiCoeffVector = zeros(1,nUniqueCoefficients);
for i = 1:nUniqueCoefficients
    uniquePhiCoeffVector(i) = prod(sumsOfColumnProducts(uniqueCoefficientIndexMatrix(i,:)));
end
phiCoeffVector = uniquePhiCoeffVector(J);
productMatrix  = phiCoeffVector'.*prod(  1+((phiVector(oneMatrix,:)-1).*binaryMatrixNumberOfPhis)  ,2);

end
%%*****************************************************************************
%%*** END OF THE FUNCTION *****************************************************
%%*****************************************************************************
%%----------------------------------------------------------------------
%% ----------------------------------------
%%----------------------------------------------------------------------








