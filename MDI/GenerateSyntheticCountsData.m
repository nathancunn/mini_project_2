%%(RSS, 17/4/13)
%%Script to generate synthetic counts data
%% - hardwired for 100 items, in 5 equal sized clusters
%%
%%----------------------------------------------------------------------
%% DEFINE USEFUL VALUES ------------------------------------------------
%%----------------------------------------------------------------------
nClusters        = 5;
nItemsPerCluster = 20;
nDataItems       = nClusters * nItemsPerCluster;
nSignalFeatures  = 20;
nNoiseFeatures   = 5;
nFeatures        = nSignalFeatures + nNoiseFeatures;
outputPath       = '~/MatlabCode/MDI/Data/';
%%----------------------------------------------------------------------
%% MEAN MATRIC FOR SYNTHETIC DATA --------------------------------------
%%----------------------------------------------------------------------
working      = ceil(nClusters*(1:nSignalFeatures)/nSignalFeatures);
signalMeans  = [];
for i=1:nClusters
  shiftVal    = i*round(nSignalFeatures/nClusters);
  for j=1:nItemsPerCluster
    signalMeans = [signalMeans; circshift(working, [0, shiftVal])];
  end
end
meanMatrix  = ones(nDataItems, nFeatures);
%%REPLACE NOISE WITH SIGNAL FOR THE RELEVANT FEATURES
meanMatrix(:, 1:nSignalFeatures) = 5*signalMeans;
imagesc(meanMatrix)
%%----------------------------------------------------------------------
%% CONSTRUCT ITEM, FEATURE NAMES ---------------------------------------
%%----------------------------------------------------------------------
itemNames    = cell(1, nDataItems);
featureNames = cell(1, nFeatures);
for i=1:nDataItems
  itemNames{i} = ['DataItem', num2str(i)];
end
for i=1:nSignalFeatures
  featureNames{i} = ['SignalFeature', num2str(i)];
end
for i=1:nNoiseFeatures
  featureNames{i+nSignalFeatures} = ['NoiseFeature', num2str(i)];
end
%%----------------------------------------------------------------------
%% WRITE DATA SETS OUT TO FILE -----------------------------------------
%%----------------------------------------------------------------------
for k=1:2
  outputFile = [outputPath, 'PoissonTestData', num2str(k),'.csv'];
  %%GENERATE DATA REALISATION
  outputData = poissrnd(meanMatrix);
  imagesc(outputData)
  %%WRITE HEADER TO FILE
  outputString = '';
  for i=1:nFeatures, outputString = [outputString, ', ', featureNames{i}, '_D', num2str(k)];, end
  outputString = [outputString, '\n'];
  fid = fopen(outputFile, 'wt');
  fprintf(fid, outputString); 
  %%WRITE EACH LINE OUT TO FILE
  for i=1:nDataItems
    outputString = itemNames{i};
    for j=1:nFeatures
      outputString = [outputString, ', ', num2str(outputData(i,j))];
    end
    outputString = [outputString, '\n'];
    fprintf(fid, outputString); 
  end
  fclose(fid);
end


%%*****************************************************************************
%%*** END OF THE SCRIPT *******************************************************
%%*****************************************************************************
%%----------------------------------------------------------------------
%% ----------------------------------------
%%----------------------------------------------------------------------

