%%Code to handle the Gaussian data type.
%%This data type treats each data point being real-valued with
%Gaussian noise.  The features are treated as independent
%measurements of the overall clustering structure.
%%
%%NOTE: The feature selection assumes that each feature is
%normalised to zero mean and unit variance
%%
%%Code originally written by Paul Kirk.
%%Additional modifications by Rich Savage.
%%
function [output1 output2] = Gaussian(input, mode, geneIndex, currentData)
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
      'empiricalMeans', [],...
      'squaredResiduals', [],...
      'mu0', [], ...
      'kappa0', [],...
      'alpha0', [],...
      'beta0', [],...
      'logZ0', [],...
      'nGenes', [], ...
      'logMarginalLikelihood', [],...
      'logicalGeneIDs', [],...
      'N', [], ...
      'featureParameters', [], ...
      'junkMeans', [], ...
      'junkVariances', []);
  %%PUT VALUES INTO THE CLUSTER STRUCTURE
  [clusterStruct.nFeatures       ] = deal(nFeatures);
  [clusterStruct.nGenesOverall   ] = deal(nGenesOverall);
  %%COMPUTE MEAN, VARIANCES FOR THE INDIFFERENCE DISTRIBUTIONS
  %%this is used for the feature selection
  [clusterStruct.junkMeans]     = deal(mean(data));
  [clusterStruct.junkVariances] = deal(std(data).^2);
  %HYPERPARAMETERS
  mu0    = 0;
  kappa0 = 0.001;
  alpha0 = 2;
  beta0  = 0.5;
  logZ0  = gammaln(alpha0) - (alpha0*log(beta0)) + (0.5*log(2*pi)) - (0.5*log(kappa0));
  %%STORE THE HYPERPARAMETER VALUES
  [clusterStruct.mu0]    = deal(mu0);
  [clusterStruct.kappa0] = deal(kappa0);
  [clusterStruct.alpha0] = deal(alpha0);
  [clusterStruct.beta0 ] = deal(beta0);
  [clusterStruct.logZ0 ] = deal(logZ0);
  %%INITIALISE VARIOUS OTHER VALUES OF INTEREST
  [clusterStruct.logMarginalLikelihood] = deal(0);
  [clusterStruct.nGenes]                = deal(0);
  [clusterStruct.logicalGeneIDs]        = deal(sparseVector);
  [clusterStruct.featureParameters]     = deal(featureParameters);
  % Initialise Context 1 clusters:
  nStartingClusters = min(ceil(log(nGenesOverall)), maxNumberOfComponents);
  % Why is this done?
  % nStartingClusters = nStartingClusters + randi(nStartingClusters);
  clusterIDs        = random('unid', nStartingClusters, 1, nGenesOverall); %row vector
  
  for i = 1:maxNumberOfComponents
    logicalIndices                   = clusterIDs == i;
    indices                          = find(logicalIndices);
    nGenesInCluster                  = length(indices);
    dataInCluster                    = sparseMatrix;
    dataInCluster(indices,:)         = data(logicalIndices,:);
    currentCluster                   = clusterStruct(i);
    currentCluster.nGenes            = nGenesInCluster;
    
    if(nGenesInCluster == 0)
      empiricalMeans               = zeros(1,nFeatures);
    elseif(nGenesInCluster == 1)
      empiricalMeans               = dataInCluster(indices,:);
    else
      empiricalMeans               = mean(dataInCluster(indices,:));
    end
    repeatedMeans                    = ones(nGenesInCluster,1)*empiricalMeans;
    
    if(nGenesInCluster == 1)
      squaredResiduals             = zeros(1,nFeatures); 
    else
      squaredResiduals             = sum((dataInCluster(indices,:) - repeatedMeans).^2);
    end
    currentCluster.empiricalMeans    = empiricalMeans;
    currentCluster.squaredResiduals  = squaredResiduals;
    currentCluster.logicalGeneIDs    = logicalIndices;
    currentCluster.N                 = nFeatures*nGenesInCluster;
    if (nGenesInCluster > 0)
      currentCluster  = Gaussian(currentCluster, 'marginal');
    end
    clusterStruct(i) = currentCluster;
  end
  
  output1 = clusterStruct;
  output2 = clusterIDs;
  
 case 'marginal'
  
  
  %%NOTE:  Could output from this option the values required for
  %feature selection
  %%i.e. compute all switched-off/on marginal likelihoods and
  %return the log-ratios
  %%
  %%Actually, this doesn't quite work :-(
  %%What would work is returning the vector of all switched-on
  %marginal likelihoods for the current cluster.
  %%Then this routine can be called for each cluster by the feature
  %selection code.
  %%This is only useful while the features are conditionally
  %independent, but that should be okay.
  
  %%we *do* want to benefit from the array of Z0 values being
  %precomputed..
  
  %%WRT switched-off features, we still need to include their
  %effect here!  How best to do this...?
  %%Maybe we need to stick with the junk probabilities idea?
  %%(approximation to make it easy to handle)
  %%
  %%This then opens up the possibility of returning a vector of
  %likelihood ratios from this routine for the feature selection :-)
  
  
  %%This computes the marginal likelihood for a single cluster
  %%EXTRACT USEFUL VALUES FROM THE INPUT 
  %%(including feature selection)
  keep             = find(input.featureParameters);
  discard          = find(input.featureParameters==0);
  nSwitchedOn      = sum(input.featureParameters);
  nSwitchedOff     = sum(input.featureParameters==0);
  empiricalMeans   = input.empiricalMeans;
  squaredResiduals = input.squaredResiduals;
  nGenes           = input.nGenes;
  mu0              = input.mu0;
  kappa0           = input.kappa0;
  alpha0           = input.alpha0;
  beta0            = input.beta0;
  logZ0            = input.logZ0;
  junkMeans        = input.junkMeans;
  junkVariances    = input.junkVariances;
  junkPrecisions   = 1./junkVariances;
  output2          = zeros(1, input.nFeatures);
  %%COMPUTE CONTRIBUTION FOR SWITCHED-ON FEATURES
  kappaN = kappa0 + nGenes;
  alphaN = alpha0 + (nGenes/2);
  betaN  = beta0  + (squaredResiduals/2) + (kappa0*nGenes*((empiricalMeans- mu0).^2))/(2*kappaN);
  logZn  = gammaln(alphaN) - (alphaN*log(betaN)) + (0.5*log(2*pi)) - (0.5*log(kappaN));
  logEv  = logZn - logZ0 - ((nGenes/2)*log(2*pi));
  %%COMPUTE SWITCHED-OFF LOG-EVIDENCE VALUES
  sumX       = nGenes           * empiricalMeans;
  sumSquares = squaredResiduals + nGenes * empiricalMeans.^2; 
  junkLogEv  = sumSquares - 2*junkMeans.*sumX + nGenes*junkMeans.^2;
  junkLogEv  = -0.5 * junkLogEv  ./ junkVariances;
  junkLogEv  = junkLogEv - 0.5*nGenes*log(2*pi*junkVariances);
  %%priors
  %%don't use these, as our switched-off model assumes that the
  %mean, variance are known quantities...
%  junkLogEv  = junkLogEv - logZ0;
%  junkLogEv  = junkLogEv + (alpha0-0.5)*log(junkPrecisions);
%  working    = kappa0*(junkMeans-mu0).^2 + 2*beta0;
%  junkLogEv  = junkLogEv - 0.5*junkPrecisions.*working; 

  logEvRatio = logEv - junkLogEv;
  output2    = logEvRatio;%%optional output for use
                                %with feature selection  
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
  %%THERE IS A FASTER WAY TO HANDLE THE MEANS, SQUARED RESIDUALS
  %%(use the profiler to see the relative speed)
  %%really, we should just propagate the sum of x, x^2
  %%but we can reverse engineer from the following:
  %% - empiricalMeans
  %% - squaredResiduals        
  %%
  %%to compute the new squaredResiduals, we use an identity:
  %%sqRes = sum(x^2) - sum(x)^2 / N
  %see eg.:
  %% http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
  %%  
  if (nGenesInCluster==0)
    newMeans      = zeros(1,output1.nFeatures);
    newSquared    = zeros(1,output1.nFeatures);
  else
    extractedData = currentData(geneIndex,:);
    oldSum        = input.nGenes           * input.empiricalMeans;
    oldSumSquares = input.squaredResiduals + oldSum.^2 / input.nGenes; 
    newSum        = oldSum                 - extractedData;
    newSumSquares = oldSumSquares          - extractedData.^2;
    newMeans      = newSum / nGenesInCluster;
    newSquared    = newSumSquares - nGenesInCluster * newMeans.^2;
  end
  %%COPY THE NEW VALUES INTO THE CORRECT VARIABLES
  empiricalMeans   = newMeans;
  squaredResiduals = newSquared;
  %%STORE USEFUL VALUES; RETURN RESULTS
  output1.empiricalMeans    = empiricalMeans;
  output1.squaredResiduals  = squaredResiduals;
  output1.logicalGeneIDs    = logicalGeneIDs;
  output1.nGenes            = nGenesInCluster;
  output1.N                 = nGenesInCluster*output1.nFeatures;
  if (nGenesInCluster == 0)
    output1.logMarginalLikelihood = 0;
  else
    output1 = Gaussian(output1, 'marginal');
  end
  
 case 'addGene'
  output1                   = input;
  initialNGenes             = output1.nGenes;
  nGenesInCluster           =  initialNGenes + 1;
  logicalGeneIDs            = output1.logicalGeneIDs;
  logicalGeneIDs(geneIndex) = true;        
  %%THERE IS A FASTER WAY TO HANDLE THE MEANS, SQUARED RESIDUALS
  %%(use the profiler to see the relative speed)
  %%really, we should just propagate the sum of x, x^2
  %%but we can reverse engineer from the following:
  %% - empiricalMeans
  %% - squaredResiduals        
  %%
  %%to compute the new squaredResiduals, we use an identity:
  %%sqRes = sum(x^2) - sum(x)^2 / N
  %see eg.:
  %% http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
  %%   
  addedData     = currentData(geneIndex,:);
  oldSum        = input.nGenes           * input.empiricalMeans;
  oldSumSquares = input.squaredResiduals + input.nGenes * input.empiricalMeans.^2; 
  newSum        = oldSum                 + addedData;
  newSumSquares = oldSumSquares          + addedData.^2;
  newMeans      = newSum / nGenesInCluster;
  newSquared    = newSumSquares - nGenesInCluster * newMeans.^2;
  %%COPY THE NEW VALUES INTO THE CORRECT VARIABLES
  empiricalMeans   = newMeans;
  squaredResiduals = newSquared;
  %%STORE USEFUL VALUES; RETURN RESULTS
  output1.empiricalMeans    = empiricalMeans;
  output1.squaredResiduals  = squaredResiduals;
  output1.logicalGeneIDs    = logicalGeneIDs;
  output1.nGenes            = nGenesInCluster;
  output1.N                 = nGenesInCluster*output1.nFeatures;
  
  if (nGenesInCluster == 0)
    output1.logMarginalLikelihood = 0;
  else
    output1 = Gaussian(output1, 'marginal');
  end
  
  %%----------------------------------------------------------------------
  %% RESAMPLE THE FEATURE PARAMETERS -------------------------------------
  %%----------------------------------------------------------------------       
  %%this code resamples all the feature parameters
  %%Because the Gaussian model assumes conditional independence
  %%between the features, we can vectorise the computations here.
  %%In essence, each feature is considered relative to the current
  %%clustering partition.
  %%Note that 'input' here needs to be the whole array of
  %ClusterStructs
  %%
  %%We compare the regular "switched-on" posterior against a
  %model where (in a given feature) all items belong to the
  %same Gaussian distribution, with known mean and variance
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
    [dum, newLogEvRatio] = Gaussian(input(currentIndex), 'marginal');
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



     
%     mu0    = input(1).mu0;
%     kappa0 = input(1).kappa0;
%     alpha0 = input(1).alpha0;
%     beta0  = input(1).beta0;
%     logZ0  = input(1).logZ0;
%     mun    = ((kappa0*mu0)+(nGenes*empiricalMeans))/(kappa0 + nGenes);
%     kappan = kappa0 + nDataItems;
%     alphan = alpha0 + (nDataItems/2);
%   betan  = beta0  + (sigmaSquared/2) + (kappa0*nDataItems*((meanValues- mu0).^2))/(2*(kappa0 + nDataItems));
%     logZn  = log(gamma(alphan)) - (alphan*log(betan)) + (0.5*log(2*pi)) - (0.5*log(kappan));

%     junkLogProbs = logZn - logZ0 - ((nDataItems/2)*log(2*pi));
     
%     pause

     
     


        %%SWITCHED-ON FEATURES FOR SINGLE NOISE VARIANCE MODEL
%        alphaN  = alpha0 + 0.5*nGenes*nSwitchedOn;
%        kappaN  = kappa0 + nGenes;
%        betaN   = beta0  + 0.5*sum(sumSquares(keep));
%        betaN   = betaN  + nGenes * kappa0 * sum(empiricalMeans(keep).^2) / (2*kappaN);
%        logLike = -0.5*nGenes*nSwitchedOn*log(2*pi);
%        logLike = logLike + 0.5*nSwitchedOn*log(kappa0/kappaN);
%        logLike = logLike + gammaln(alphaN);
%        logLike = logLike + alpha0*log(beta0);
%        logLike = logLike - gammaln(alpha0);
%        logLike = logLike - alphaN*log(betaN);        






