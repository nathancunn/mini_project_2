alldata = importdata('Data/GaussianTestData1.csv', ',',1);
alldata = importdata('Data/GaussianTestData1.csv', ',',1);
dataForThisContext = alldata.data;
numbofparts = 100;
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
D = 2;
%xi = ones(d, numbofparts);
s = zeros(nGenes, numbofparts, D);
prob = ones(1, D) / D;


proposedClustersForThisContext = clusterContainer(1).clusterStruct
for ind = 1:nGenes
    proposedClustersForThisContext = AddRemoveItem('removeGene', proposedClustersForThisContext, s(ind, 1), ind, dataForThisContext, 1);
end
clusterContainer(1).clusterStruct = proposedClustersForThisContext
particle_clusters = cell(1, numbofparts);
particle_clusters(1, :) = {clusterContainer}
clusterContainer(2:200) = clusterContainer(1);
proposedClusterContainer = clusterContainer;
part_inds = ((0:(nGenes - 1)) * numberOfComponents) + 1
s = zeros(nGenes, numbofparts, D);
prob = ones(1, D) / D;
logweight = zeros(D, numbofparts);
for i = 1:nGenes
    disp(['i = ' num2str(i)])
    for m = 1:numbofparts
        for d = 1:D;
            
            if(i == 1)
                proposedClusterContainer(D * (m - 1) + d).clusterStruct = AddRemoveItem('addGene', proposedClusterContainer(D * (m - 1) + d).clusterStruct, 1, i, dataForThisContext, 1);
            else
                logprob = zeros(1, 2);
                for ind = 1:2
                    prop = AddRemoveItem('addGene', proposedClusterContainer(D * (m - 1) + d).clusterStruct, ind, i, dataForThisContext, 1);
                    logprob(1,ind) = sum([exp(prop.logMarginalLikelihood]));
                end
                fprob = cumsum(prob .* exp(logprob - max(logprob)));
                fprob = fprob/fprob(end);
                u1 = rand;
                sstar = 1;
                while ( fprob(sstar) < u1 )
                    sstar = sstar + 1;
                end
                logweight(d, m) = logweight(d, m) + log(fprob(end)) + max(logprob);
                proposedClusterContainer(D * (m - 1) + d).clusterStruct = AddRemoveItem('addGene', proposedClusterContainer(D * (m - 1) + d).clusterStruct, sstar, i, dataForThisContext, 1);

            end
            s(i,m,d) = sstar;
        end
        
    end
    clusterContainer = proposedClusterContainer(1:D);
end

dataForThisContext = clusterContainer(1).data;
currentCluster = clusterContainer(1).clusterStruct;
Gaussian(currentCluster(2), 'marginal')
Gaussian(s(:, 1), 'marginal')
Gaussian(test, 'marginal', 1, dataForThisContext)

maxNumberOfComponents = 2;
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
  
  nGenesOverall         = size(data,1);
  nFeatures             = size(data,2);
  maxNumberOfComponents = 2;
  sparseVector          = false(1,nGenesOverall);
  sparseMatrix          = zeros(nGenesOverall,nFeatures);
  featureParameters     = input.featureParameters;
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
  
  clusterContainer(1).clusterStruct = Gaussian(clusterStruct, 'init')
  Gaussian(clusterContainer(1), 'addGene', 1);
  
  proposedClustersForThisContext = clusterContainer(1).clusterStruct;
for ind = 1:nGenes
    proposedClustersForThisContext = AddRemoveItem('removeGene', proposedClustersForThisContext, s(ind, 1), ind, dataForThisContext, 1);
end
proposedClustersForThisContext = AddRemoveItem('removeGene', proposedClustersForThisContext, 1, 1, dataForThisContext, 1);
  proposedClustersForThisContext = AddRemoveItem('addGene', proposedClustersForThisContext, 1, 1, dataForThisContext, 1);
  proposedClustersForThisContext = AddRemoveItem('addGene', proposedClustersForThisContext, 1, 1, dataForThisContext, 1);
