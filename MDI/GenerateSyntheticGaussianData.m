%%Script to generate some synthetic Gaussian data with which to
%test MDI
%%----------------------------------------------------------------------
%% DEFINE USEFUL VALUES ------------------------------------------------
%%----------------------------------------------------------------------
nDataItems      = 200;
nClusters       = 3;
nSignalFeatures = 100;
nNoiseFeatures  = 50;
signalLevel     = 30;
nFeatures       = nSignalFeatures + nNoiseFeatures;
outputPath      = 'Data/synth_data/';
outputFile1     = [outputPath, 'aGaussianTestData1.csv'];
outputFile2     = [outputPath, 'aGaussianTestData2.csv'];
%%----------------------------------------------------------------------
%% CONSTRUCT SYNTHETIC DATA --------------------------------------------
%%----------------------------------------------------------------------
data        = zeros(nDataItems, nFeatures);
signalArray = ceil((1:nDataItems)*nClusters/nDataItems);
signalArray = signalArray - mean(signalArray);
signalArray = signalLevel * signalArray / nClusters;
%%RANDOMLY SHUFFLE THE ITEMS
[dum, sortIndex] = sort(rand(1, nDataItems));
signalArray      = signalArray(sortIndex);
%%ASSIGN THE SIGNAL FEATURES
for i=1:nSignalFeatures
  data(:, i) = data(:, i) + signalArray';
end
imagesc(data)
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
  outputFile = [outputPath, 'GaussianTestData', num2str(k),'.csv'];
  %%ADD A NOISE REALISATION TO THE DATA
  outputData = data + randn(nDataItems, nFeatures);  
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

