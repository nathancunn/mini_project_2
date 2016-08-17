%%(6/1/12 Rich Savage)
%%Function to find the consensus clustering structure from an MDI run.
%%The MCMC files are structured as follows:
%%
%%> 2) a CSV file that contains the sampled values. Let K=nDatasets. The first K 
%%> columns provide the sampled values for the K mass parameters associated with 
%%> each Dirichlet prior. The next K(K-1)/2 columns provide the "phi" (dataset 
%%> association) parameters. These are in the order phi_{1,2}, phi_{1,3}, ..., 
%%> phi_{1,K}, phi_{2,3}, phi_{2, 4}, ... and so on. The next p columns provide 
%%> the component allocation labels for the p genes in Dataset 1. The next p 
%%> columns contain the component allocation labels for the p genes in Dataset 2. 
%> ... and so on.
%%
function outputIDs = GenerateClusteringPartition(outputPath, mcmcFiles, nDataTypes, burnInFraction, typeIndex)
%%----------------------------------------------------------------------
%% HANDLE THE OPTIONAL typeIndex INPUT ---------------------------------
%%----------------------------------------------------------------------
if exist('typeIndex', 'var')==0
  typeIndex = 1:nDataTypes;
end
%%----------------------------------------------------------------------
%% FIND THE ITEM NAMES -------------------------------------------------
%%----------------------------------------------------------------------
[dum, paramNames] = ReadInMcmcFile(mcmcFiles{1});
working           = regexp(paramNames, 'Dataset1_');
keep              = find(cellfun('length', working));
itemNames         = paramNames(keep);
nDataItems        = length(itemNames);
nActualTypes      = length(typeIndex);
for i=1:nDataItems
  itemNames{i} = itemNames{i}(10:end);
end
%%----------------------------------------------------------------------
%% FIND USEFUL VALUES --------------------------------------------------
%%----------------------------------------------------------------------
nParams           = nDataTypes + nDataTypes*(nDataTypes-1)/2;  %%the subset of parameters to display
nMcmcFiles        = length(mcmcFiles);
plotFile          = [outputPath, '_nClusterDiagnostics.pdf'];
outputFile        = [outputPath, '_clusteringPartition.csv'];
similarityFile    = [outputPath, '_similarityMatrix.jpg'];
sampleCounter     = 0;
similarityArray   = zeros(nDataItems, nDataItems, nDataTypes);
nClusterArray     = [];
nFusionPairwise   = nDataTypes*(nDataTypes-1)/2; 
fusionStates      = zeros(nDataItems, nFusionPairwise);
%%----------------------------------------------------------------------
%% FIND THE POSTERIOR SIMILARITY MATRIX;  PLOT nCluster TIME SERIES ----
%%----------------------------------------------------------------------
figure(1)
subplot(nDataTypes, 2, 1)
for i=1:nMcmcFiles
  %%READ IN THE CURRENT DATA FILE
  data            = ReadInMcmcFile(mcmcFiles{i});
  nSamples        = size(data,1);
  keep            = ceil(burnInFraction*nSamples):nSamples;
  data            = data(keep, :);
  nSamples        = length(keep);
  nClusterCurrent = zeros(nSamples, nDataTypes);
  %%FOR EACH PAIR OF DATA TYPES, FIND THE FUSION STATES
  fusionCounter     = 0;
  for j=1:(nDataTypes-1)
    for k=(j+1):nDataTypes
      fusionCounter = fusionCounter+1;
      index_j       = (1:nDataItems) + nParams + (j-1)*nDataItems; 
      index_k       = (1:nDataItems) + nParams + (k-1)*nDataItems; 
      data_j        = data(:, index_j);
      data_k        = data(:, index_k);
      fusionStates(:, fusionCounter) = fusionStates(:, fusionCounter) + sum(data_j==data_k)';
    end
  end
  %%FOR EACH DATA TYPE, ADD TO THE POSTERIOR SIMILARITY
  for j=1:nDataTypes
    %%FIND THE CLUSTER IDS FOR THE CURRENT DATA TYPE
    index       = 1:nDataItems;
    index       = index + nParams + (j-1)*nDataItems; 
    currentData = data(:, index);
    %%CONTRIBUTION TO THE POSTERIOR SIMILARITY MATRIX
    for k=1:nDataItems
      currentIDs              = currentData(:,k);
      currentIDs              = currentIDs * ones(1, nDataItems);
      idMatch                 = currentIDs==currentData;
      similarityArray(k,:,j) = similarityArray(k,:,j) + sum(idMatch, 1);
    end    
    %%COUNT THE NUMBER OF CLUSTERS 
    for k=1:nSamples
      nClusterCurrent(k,j) = length(unique(currentData(k,:)));
    end
    %%PLOT THE nCluster TIME SERIES
    if ismember(j, typeIndex)
      subplot(nActualTypes, 2, find(typeIndex==j))
      plot(nClusterCurrent(:,j))
      hold all
    end
  end
  %%UPDATE SAMPLE COUNTER; STORE THE nCluster VALUES
  sampleCounter = sampleCounter + nSamples;
  nClusterArray = [nClusterArray; nClusterCurrent];
end
%%----------------------------------------------------------------------
%% SAVE FIGURE TO FILE -------------------------------------------------
%%----------------------------------------------------------------------
maxIndex = max(nClusterArray(:));
for i=1:nDataTypes
  if ismember(i, typeIndex)
    currentIndex = find(typeIndex==i);
    subplot(nActualTypes, 2, currentIndex)
    title(['dataType ', num2str(i)], 'FontSize', 6, 'Interpreter', 'none')
    xlabel('sample')
    ylabel('nClusters')
    ylim([0, 1.1*maxIndex])
    hold off
    subplot(nActualTypes, 2, currentIndex+nActualTypes)
    title(['dataType ', num2str(i)], 'FontSize', 6, 'Interpreter', 'none')
    hist(nClusterArray(:,i), 1:maxIndex)
    hold off
  end 
end
saveas(gcf, plotFile, 'pdf')
subplot(1,1,1)
%%----------------------------------------------------------------------
%% OPTION TO MAKE A CUT ON P(FUSION) -----------------------------------
%%----------------------------------------------------------------------
fusionStates    = fusionStates ./ sampleCounter;
meanFusionProbs = mean(fusionStates, 2);
%keep            = find(meanFusionProbs>=median(meanFusionProbs));         
%similarityArray = similarityArray(keep, keep, :);
%fusionStates    = fusionStates(keep);
%itemNames       = itemNames(keep);
%nDataItems      = length(itemNames);
%%----------------------------------------------------------------------
%% COMPUTE SIMILARITY MATRIX (OPTIONALLY FROM A SUBSET OF DATA TYPES) --
%%----------------------------------------------------------------------
similarityArray  = similarityArray                        ./ sampleCounter;
similarityMatrix = sum(similarityArray(:,:,typeIndex), 3) ./ length(typeIndex);
working          = nClusterArray(:, typeIndex);
nMapClusters     = round(mean(working(:)));
%%----------------------------------------------------------------------
%% EXTRACT THE CLUSTERING PARTITION ------------------------------------
%%----------------------------------------------------------------------
%%Medvedovic used hierarchical clustering, using the posterior
%similarities as distances.  Do this, taking nClusters
%from the posterior distribution
%%
%%find a single, overall similarity matrix
%similarityArray  = similarityArray         ./ sampleCounter;
%similarityMatrix = sum(similarityArray, 3) ./ nDataTypes;
%nMapClusters     = round(mean(nClusterArray(:)));
%%FORMAT THE DISTANCES FOR MATLAB LIBRARY FUNCTIONS
%%'distance' is defined as 1-P(similarity)
distances = [];
for i=1:(nDataItems-1)
  distances = [distances, 1 - similarityMatrix(i, (i+1):end)];
end
%%FIND CLUSTERING PARTITION 
working    = linkage(distances, 'complete');
clusterIDs = cluster(working, 'maxclust', nMapClusters);
outputIDs  = clusterIDs;%%preserve the original order for the output
%%SORT BY THE CLUSTERING PARTITION
[clusterIDs, sortIndex] = sort(clusterIDs);
similarityMatrix        = similarityMatrix(sortIndex, sortIndex);
fusionStates            = fusionStates(sortIndex, :);
itemNames               = itemNames(sortIndex);
%%----------------------------------------------------------------------
%% WRITE THE CLUSTERING PARTITION TO FILE ------------------------------
%%----------------------------------------------------------------------
fid = fopen(outputFile, 'wt');
for i=1:nDataItems
  fprintf(fid, [itemNames{i}, ',', num2str(clusterIDs(i)), '\n']); 
end
fclose(fid);
%%----------------------------------------------------------------------
%% PLOT THE SIMILARITY MATRIX TO FILE ----------------------------------
%%----------------------------------------------------------------------
%%SHOW A MEASURE OF THE FUSION STATE
%%try the sum over pairwise fusion states
%%(this may not be sensible for nDataTypes>2)
meanFusionStates = mean(fusionStates, 2);
fusionLabels     = repmat({''}, 1, nDataItems);
nIncrements      = 10;
for i=1:nIncrements
  index               = find(meanFusionStates>(i/nIncrements));
  fusionLabels(index) = strcat(fusionLabels(index), '-');
end
%%PLOT IMAGES
imagesc(similarityMatrix)
colorbar
index                  = find(mod(clusterIDs, 2));
partitionLabels        = repmat({'.'}, 1, nDataItems);
partitionLabels(index) = {'l'};
set(gca, 'XTick',      1:nDataItems)
set(gca, 'XTickLabel', partitionLabels)
xlabel('(partition)')
set(gca, 'YTick',      1:nDataItems)
set(gca, 'YTickLabel', fusionLabels)
ylabel('(fusion state)')
saveas(gcf, similarityFile, 'jpg')
%%ALSO PLOT FOR EACH DATA TYPE
if length(typeIndex)>1
  for i=1:length(typeIndex)
    imagesc(similarityArray(sortIndex,sortIndex,typeIndex(i)))
    colorbar
    set(gca, 'XTick',      1:nDataItems)
    set(gca, 'XTickLabel', partitionLabels)
    xlabel('(partition)')
    set(gca, 'YTick',      1:nDataItems)
    set(gca, 'YTickLabel', fusionLabels)
    ylabel('(fusion state)')
    saveas(gcf, [outputPath, '_similarityArray', num2str(typeIndex(i)), '.jpg'], 'jpg')
  end
end
%%----------------------------------------------------------------------
%% WRITE THE SIMILARITY MATRIX TO A CSV FILE ---------------------------
%%----------------------------------------------------------------------
%%can get the itemNames from the partition file
csvwrite([outputPath, '_similarityMatrix.csv'], similarityMatrix)
%%*****************************************************************************
%%*** END OF THE FUNCTION *****************************************************
%%*****************************************************************************
%%----------------------------------------------------------------------
%% ----------------------------------------
%%----------------------------------------------------------------------



