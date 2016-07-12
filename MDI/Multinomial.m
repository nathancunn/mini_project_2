%%Code to handle the multinomial data type.
%%This data type treats each feature as being drawn from a
%multinomial distribution, and forms the overal marginal likelihood
%by taking the product of the multinomials for each feature
%%
%%Code originally written by Paul Kirk.
%%Additional modifications by Rich Savage.
%%
function [output1 output2] = Multinomial(input, mode, geneIndex, dataForCurrentGene)
switch mode
      case 'init'
        %%EXTRACT USEFUL INPUT VALUES
        data                  = input.data;
        nGenesOverall         = input.nGenes;
        nFeatures             = input.nFeatures;
        maxNumberOfComponents = input.maxNumberOfComponents;
        sparseVector          = false(1,nGenesOverall);
        sparseMatrix          = zeros(nGenesOverall,nFeatures);
        featureParameters     = input.featureParameters;
        %%GENERATE SOME PLOTS OF THE DATA
        %%this is for diagnostic purposes
        %%the user should probably always at least glance at the
        %input data!
        figure(1), imagesc(data), colorbar
        figure(2), hist(data(:))
        pause(3)
        %%DEFINE THE CLUSTER STRUCTURE
        clusterStruct(1,maxNumberOfComponents) = struct(...
            'nFeatures', [], ...
            'nGenesOverall', [], ...
            'dataLevels', [],...
            'dataCountIndexHelper', [],...
            'beta', [], ...
            'sumBeta',[],...
            'nGenes', [], ...
            'logMarginalLikelihood', [],...
            'dataCounts', [], ...
            'logicalGeneIDs', [], ...
            'N', [], ...
            'featureParameters', [], ...
            'junkLogProbs', []);
          %%FIND THE JUNK PROBABILITIES
          %%these should represent the distribution of data values for each feature
          %%hence, this should be a (nDataValues*nFeatures) array
          nDataValues = length(unique(data));
          junkProbs   = zeros(nDataValues, nFeatures);
          for i=1:nFeatures
            currentProbs    = histc(data(:,i), 1:nDataValues,1);
            currentProbs    = max(currentProbs, 1);%%avoid zeros here
            currentProbs    = currentProbs / sum(currentProbs);
            junkProbs(:, i) = currentProbs;
          end
          [clusterStruct.junkLogProbs] = deal(log(junkProbs));
        %%PUT VALUES INTO THE CLUSTER STRUCTURE
        [clusterStruct.nFeatures       ]      = deal(nFeatures);
        [clusterStruct.nGenesOverall   ]      = deal(nGenesOverall);
        dataLevels                            = unique(data);
        dataLevels                            = reshape(dataLevels, 1, length(dataLevels)); %ensure we have a row vector
        [clusterStruct.dataLevels]            = deal(dataLevels);
        nLevels                               = length(dataLevels);
        dataCountIndexHelper                  = 0:nLevels:(nLevels*(nFeatures-1));
        [clusterStruct.dataCountIndexHelper ] = deal(dataCountIndexHelper);
        dataProportions                       = histc(data(:), dataLevels);        
        dataProportions                       = dataProportions/sum(dataProportions);
        hyperParameters                       = length(dataProportions)*0.5*dataProportions;
        beta                                  = repmat(hyperParameters, 1, nFeatures);
        [clusterStruct.beta]                  = deal(beta);
        [clusterStruct.sumBeta]               = deal(sum(beta, 1));
        [clusterStruct.logMarginalLikelihood] = deal(0);
        [clusterStruct.nGenes]                = deal(0);
        [clusterStruct.logicalGeneIDs]        = deal(sparseVector);
        [clusterStruct.featureParameters]     = deal(featureParameters);
        % Initialise Context 1 clusters:
        nStartingClusters = min(ceil(log(nGenesOverall)), maxNumberOfComponents);
        clusterIDs        = random('unid', nStartingClusters, 1, nGenesOverall); 
        
        for i = 1:maxNumberOfComponents
            logicalIndices                   = clusterIDs == i;
            indices                          = find(logicalIndices);
            nGenesInCluster                  = length(indices);
            dataInCluster                    = sparseMatrix;
            dataInCluster(indices,:)         = data(logicalIndices,:);
            currentCluster                   = clusterStruct(i);
            currentCluster.nGenes            = nGenesInCluster;
            currentCluster.dataCounts        = histc(dataInCluster(indices,:),dataLevels);
            currentCluster.logicalGeneIDs    = logicalIndices;
            currentCluster.N                 = nFeatures*nGenesInCluster;
            if (nGenesInCluster > 0)
                currentCluster  = Multinomial(currentCluster, 'marginal'); 
            end
            clusterStruct(i) = currentCluster;
        end
        
        output1 = clusterStruct;
        output2 = clusterIDs;
        
 case 'marginal'
  %%This computes the marginal likelihood for a single cluster
  %%EXTRACT USEFUL VALUES FROM THE INPUT 
  %%(including feature selection)
  keep         = find(input.featureParameters);
  discard      = find(input.featureParameters==0);
  nSwitchedOn  = sum(input.featureParameters);
  nSwitchedOff = sum(input.featureParameters==0);
  dataCounts   = input.dataCounts(:, keep);
  beta         = input.beta(:, keep);
  sumBeta      = input.sumBeta(keep);
  nGenes       = input.nGenes;
  %%COMPUTE THE NEW LOG MARGINAL LIKELIHOOD
  logMarginalLikelihood = sum(sum(gammaln(dataCounts+beta)));
  logMarginalLikelihood = logMarginalLikelihood + sum(gammaln(sumBeta));
  logMarginalLikelihood = logMarginalLikelihood - sum(sum(gammaln(beta)));
  logMarginalLikelihood = logMarginalLikelihood - sum(gammaln(nGenes + sumBeta));
  %%IF REQUIRED, UPDATE SWITCHED-OFF LOG-EVIDENCE VALUES
  if nSwitchedOff>0
    dataCounts_off        = input.dataCounts(:, discard);
    beta_off              = input.beta(:, discard);
    junkLogProbs          = input.junkLogProbs(:, discard);  
    newLogEv              = sum(sum((dataCounts_off+beta_off).*junkLogProbs));
    logMarginalLikelihood = logMarginalLikelihood + newLogEv;
  end
  %%STORE AND RETURN THE NEW LOG MARGINAL LIKELIHOOD
  input.logMarginalLikelihood = logMarginalLikelihood;
  output1 = input;
  
  %%CAN THIS BE DELETED???
    case 'initialiseAuxiliary'
        output1 = input;
        nGenesInCluster    = 1;
        output1.nGenes     = nGenesInCluster;
        output1.N          = input.nFeatures;
        output1.dataCounts = histc(dataForCurrentGene, input.dataLevels,1);
        output1.logicalGeneIDs = true;
        output1 = Multinomial(output1, 'marginal');
    case 'removeGene'
        output1 = input;
        %%%%%
        dataCountIndexHelper = output1.dataCountIndexHelper;
        dataCounts = output1.dataCounts;
        indices = dataCountIndexHelper+dataForCurrentGene;
        nGenesInCluster   = output1.nGenes - 1;
        dataCounts(indices) = dataCounts(indices) - 1;
        logicalGeneIDs = output1.logicalGeneIDs;
        logicalGeneIDs(geneIndex) = false;
        output1.logicalGeneIDs = logicalGeneIDs;
        output1.nGenes            = nGenesInCluster;
        output1.dataCounts        = dataCounts;
        output1.N                 = nGenesInCluster*output1.nFeatures;
        if (nGenesInCluster == 0)
            output1.logMarginalLikelihood = 0;
        else
            output1 = Multinomial(output1, 'marginal');
        end
    case 'addGene'
        output1 = input;
        initialNGenes = output1.nGenes;
        dataCountIndexHelper = output1.dataCountIndexHelper;
        
        if(initialNGenes>0)
            dataCounts = output1.dataCounts;
        else
            dataCounts = zeros(dataCountIndexHelper(2),output1.nFeatures);
        end
        indices = dataCountIndexHelper+dataForCurrentGene;
        nGenesInCluster   =  initialNGenes + 1;
        dataCounts(indices) = dataCounts(indices) + 1;
        logicalGeneIDs = output1.logicalGeneIDs;
        logicalGeneIDs(geneIndex) = true;
        output1.logicalGeneIDs = logicalGeneIDs;
        
        output1.nGenes            = nGenesInCluster;
        output1.dataCounts        = dataCounts;
        output1.N                 = nGenesInCluster*output1.nFeatures;
        if (nGenesInCluster == 0)
            output1.logMarginalLikelihood = 0;
        else
            output1 = Multinomial(output1, 'marginal');
        end
  
 %%----------------------------------------------------------------------
 %% RESAMPLE THE FEATURE PARAMETERS -------------------------------------
 %%----------------------------------------------------------------------       
 %%this code resamples all the feature parameters
 %%Because the multinomil model assumes conditional independence
 %%between the features, we can vectorise the computations here.
 %%In essence, each feature is considered relative to the current
 %%clustering partition.
 %%Note that 'input' here needs to be the whole array of
 %ClusterStructs
 %%
 %%NOTE:  feature selection here is super-fast!
 %This suggests that we might be able to include lots of features
 %%esp if we initialise the clusterIDs in a sensible way
 case 'featureSelection'
  %%FIND USEFUL VALUES
  clusterIDs   = geneIndex; %%odd name due to input convention!
  uniqueIDs    = unique(clusterIDs);
  nClusters    = length(uniqueIDs);
  nFeatures    = input(1).nFeatures;
  nDataLevels  = length(input(1).dataLevels);
  pSwitchedOff = zeros(1, nFeatures);
  logProbArray = zeros(1, nFeatures);
  junkLogProbs = input(1).junkLogProbs;
  %%FIND THE LOG PROBABILITY RATIOS
  %%loop over the occupied clusters
  %%vectorised over the features 
  for i=1:nClusters
    currentIndex = uniqueIDs(i);
    dataCounts   = input(currentIndex).dataCounts;
    sumBeta      = input(currentIndex).sumBeta;
    beta         = input(currentIndex).beta;
    nGenes       = input(currentIndex).nGenes;
    newLogEv     = sum(gammaln(dataCounts+beta));
    newLogEv     = newLogEv   + gammaln(sumBeta);
    newLogEv     = newLogEv   - sum(gammaln(beta));
    newLogEv     = newLogEv   - gammaln(nGenes + sumBeta);
    newLogEv     = newLogEv   - sum((dataCounts+beta).*junkLogProbs);%%W=0 contribution
    logProbArray = logProbArray + nGenes .* newLogEv;
  end
  %%FIND THE CONDITIONAL DISTRIBUTION FOR EACH FEATURE
  probArray    = exp(logProbArray);
  pSwitchedOff = 1 ./ (probArray+1);
  pSwitchedOn  = 1 - pSwitchedOff;
  %%GIBBS STEPS FOR EACH FEATURE
  %%this is vectorised as these steps are independent of one another
  newFeatureParameters      = unifrnd(0, 1, 1, nFeatures) > pSwitchedOff;
  [input.featureParameters] = deal(newFeatureParameters);
  output1                   = input;
end

end
%%*****************************************************************************
%%*** END OF THE FUNCTION *****************************************************
%%*****************************************************************************
%%----------------------------------------------------------------------
%% ----------------------------------------
%%----------------------------------------------------------------------









