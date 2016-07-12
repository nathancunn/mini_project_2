%%Code to handle the Poisson data type.
%%This data type treats each data point being a non-negative
%integer (count), drawn from a Poisson distribution.
%%The mean of the Poisson distribution is unknown and is
%marginalised.
%%  The features are treated as independent
%measurements of the overall clustering structure.
%%
%%Code written by Rich Savage.
%%
function [output1 output2] = Poisson(input, mode, geneIndex, currentData)
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
  figure(2), hist(data(:), 1000)
  pause(3)
  % DEFINE THE CLUSTER STRUCTURE
  clusterStruct(1,maxNumberOfComponents) = struct(...
      'nFeatures', [], ...
      'nGenesOverall', [], ...
      'sumCounts', [],...
      'sumLogGamma', [],...
      'a', [], ...
      'b', [],...
      'nGenes', [], ...
      'logMarginalLikelihood', [],...
      'logicalGeneIDs', [],...
      'featureParameters', [], ...
      'junkMeans', []);
  %%PUT VALUES INTO THE CLUSTER STRUCTURE
  [clusterStruct.nFeatures       ] = deal(nFeatures);
  [clusterStruct.nGenesOverall   ] = deal(nGenesOverall);
  %%COMPUTE MEAN, VARIANCES FOR THE INDIFFERENCE DISTRIBUTIONS
  %%this is used for the feature selection
  [clusterStruct.junkMeans]     = deal(mean(data));
  %HYPERPARAMETERS
  [clusterStruct.a] = deal(5);
  [clusterStruct.b] = deal(1);
  %%INITIALISE VARIOUS OTHER VALUES OF INTEREST
  [clusterStruct.logMarginalLikelihood] = deal(0);
  [clusterStruct.nGenes]                = deal(0);
  [clusterStruct.logicalGeneIDs]        = deal(sparseVector);
  [clusterStruct.featureParameters]     = deal(featureParameters);
  % Initialise Context 1 clusters:
  nStartingClusters = min(ceil(log(nGenesOverall)), maxNumberOfComponents);
  nStartingClusters = nStartingClusters + randi(nStartingClusters);
  clusterIDs        = random('unid', nStartingClusters, 1, nGenesOverall); %row vector
  %%INITIALISE EACH CLUSTER STRUCTURE
  for i = 1:maxNumberOfComponents
    logicalIndices                   = clusterIDs == i;
    indices                          = find(logicalIndices);
    nGenesInCluster                  = length(indices);
    dataInCluster                    = sparseMatrix;
    dataInCluster(indices,:)         = data(logicalIndices,:);
    currentCluster                   = clusterStruct(i);
    currentCluster.nGenes            = nGenesInCluster;    
    currentCluster.sumCounts         = sum(dataInCluster(indices, :), 1);
    currentCluster.sumLogGamma       = sum(gammaln(1+dataInCluster(indices,:)), 1);
    currentCluster.logicalGeneIDs    = logicalIndices;
    if (nGenesInCluster > 0)
      currentCluster  = Poisson(currentCluster, 'marginal');
    end
    clusterStruct(i) = currentCluster;
  end
  output1 = clusterStruct;
  output2 = clusterIDs;

  
 case 'marginal'
  %%This computes the marginal likelihood for a single cluster
  %%EXTRACT USEFUL VALUES FROM THE INPUT 
  %%(including feature selection)
  keep             = find(input.featureParameters);
  discard          = find(input.featureParameters==0);
  nSwitchedOn      = sum(input.featureParameters);
  nSwitchedOff     = sum(input.featureParameters==0);
  sumCounts        = input.sumCounts;
  sumLogGamma      = input.sumLogGamma;
  nGenes           = input.nGenes;
  a                = input.a;
  b                = input.b;
  junkMeans        = input.junkMeans;
  output2          = zeros(1, input.nFeatures);
  %%COMPUTE CONTRIBUTION FOR SWITCHED-ON FEATURES
  logEv = gammaln(sumCounts + a);
  logEv = logEv - (sumCounts+a)*log(nGenes+b);
  logEv = logEv - sumLogGamma;
  logEv = logEv - gammaln(a);
  logEv = logEv + a*log(b); 
  %%COMPUTE SWITCHED-OFF LOG-EVIDENCE VALUES
  junkLogEv = sumCounts .* log(junkMeans);
  junkLogEv = junkLogEv -  nGenes*junkMeans;
  junkLogEv = junkLogEv -  sumLogGamma;
  logEvRatio = logEv - junkLogEv;
  output2    = logEvRatio;%%optional output for use with feature selection  
  if nSwitchedOff>0
    logEv(discard) = junkLogEv(discard);
  end
  %%COMPUTE, STORE AND RETURN THE NEW LOG MARGINAL LIKELIHOOD; 
  input.logMarginalLikelihood = sum(logEv);
  output1                     = input;
  
 case 'removeGene'
  output1                   = input;
  nGenesInCluster           = output1.nGenes - 1;
  logicalGeneIDs            = output1.logicalGeneIDs;
  logicalGeneIDs(geneIndex) = false;
  %%COMPUTE UPDATED STATISTICS
  if (nGenesInCluster==0)
    newSumCounts   = zeros(1,output1.nFeatures);
    newSumLogGamma = zeros(1,output1.nFeatures);
  else
    extractedData  = currentData(geneIndex,:);
    newSumCounts   = input.sumCounts   - extractedData;
    newSumLogGamma = input.sumLogGamma - gammaln(1+extractedData);
  end
  %%STORE USEFUL VALUES; RETURN RESULTS
  output1.sumCounts        = newSumCounts;
  output1.sumLogGamma    = newSumLogGamma;
  output1.logicalGeneIDs = logicalGeneIDs;
  output1.nGenes         = nGenesInCluster;
  if (nGenesInCluster == 0)
    output1.logMarginalLikelihood = 0;
  else
    output1 = Poisson(output1, 'marginal');
  end
  
 case 'addGene'
  output1                   = input;
  initialNGenes             = output1.nGenes;
  nGenesInCluster           =  initialNGenes + 1;
  logicalGeneIDs            = output1.logicalGeneIDs;
  logicalGeneIDs(geneIndex) = true;        
  %%COMPUTE UPDATED STATISTICS
  addedData      = currentData(geneIndex,:);
  newSumCounts   = input.sumCounts   + addedData;
  newSumLogGamma = input.sumLogGamma + gammaln(1+addedData);
  %%STORE USEFUL VALUES; RETURN RESULTS
  output1.sumCounts      = newSumCounts;
  output1.sumLogGamma    = newSumLogGamma;
  output1.logicalGeneIDs = logicalGeneIDs;
  output1.nGenes         = nGenesInCluster;
  if (nGenesInCluster == 0)
    output1.logMarginalLikelihood = 0;
  else
    output1 = Poisson(output1, 'marginal');
  end
  
  %%----------------------------------------------------------------------
  %% RESAMPLE THE FEATURE PARAMETERS -------------------------------------
  %%----------------------------------------------------------------------       
  %%this code resamples all the feature parameters
  %%Because the model assumes conditional independence
  %%between the features, we can vectorise the computations here.
  %%In essence, each feature is considered relative to the current
  %%clustering partition.
  %%Note that 'input' here needs to be the whole array of
  %ClusterStructs
  %%
  %%We compare the regular "switched-on" posterior against a
  %model where (in a given feature) all items belong to the
  %same distribution, with known parameters
  %(estimated from all data in a given feature).  This is an
  %approximation to the more proper model of marginalising over an
  %unknown mean, variance; we do this to make the code more
  %tractable.  In practice, the effect is to slightly favour the
  %switched-off model a priori.
 case 'featureSelection'
  %%FIND USEFUL VALUES
  clusterIDs   = geneIndex; %%odd name due to input convention!
  uniqueIDs    = unique(clusterIDs);
  nClusters    = length(uniqueIDs);
  nFeatures    = input(1).nFeatures;
  logEvRatio   = zeros(1, nFeatures);
  %%FIND THE LOG PROBABILITY RATIOS
  %%loop over the occupied clusters
  %%vectorised over the features 
  for i=1:nClusters
    currentIndex         = uniqueIDs(i);
    [dum, newLogEvRatio] = Poisson(input(currentIndex), 'marginal');
    logEvRatio           = logEvRatio + newLogEvRatio;
  end
  %%FIND THE CONDITIONAL DISTRIBUTION FOR EACH FEATURE
  probArray    = exp(logEvRatio);
  pSwitchedOff = 1 ./ (probArray+1);
  pSwitchedOn  = 1 - pSwitchedOff;
  %%GIBBS STEPS FOR EACH FEATURE
  %%this is vectorised as these steps are independent of one another
  %%apply a prior that at least one feature must be switched on
  newFeatureParameters      = unifrnd(0, 1, 1, nFeatures) > pSwitchedOff;
  output1                   = input;
  if sum(newFeatureParameters)>0
    [output1.featureParameters] = deal(newFeatureParameters);
  end
end

end
%%*****************************************************************************
%%*** END OF THE FUNCTION *****************************************************
%%*****************************************************************************
%%----------------------------------------------------------------------
%% ----------------------------------------
%%----------------------------------------------------------------------





