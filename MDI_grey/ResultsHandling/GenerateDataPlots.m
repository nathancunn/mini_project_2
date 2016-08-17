%%(6/1/12 Rich Savage)
%%Function to generate data plots, sorted by clustering partition
%and (optionally) feature selection.


%%do we want all plots on a single page?? (and hence in a single
%file)
%%this would simplify the handling of multiple data type combination


%%
function GenerateDataPlots(outputPath, dataFiles, clusterIDs, featureProbs, featureNames)
%%----------------------------------------------------------------------
%% FIND USEFUL VALUES --------------------------------------------------
%%----------------------------------------------------------------------
nDataTypes              = length(dataFiles);
nDataItems              = length(clusterIDs);
[clusterIDs, sortIndex] = sort(clusterIDs);
dimensionSize           = ceil(sqrt(nDataTypes));
outputFile              = [outputPath, '_SortedData.pdf'];
%%----------------------------------------------------------------------
%% REMOVE PREFIX FROM THE FEATURE NAMES --------------------------------
%%----------------------------------------------------------------------
for i=1:length(featureNames)
  [dum remains]   = strtok(featureNames{i}, '_');
  featureNames{i} = remains(2:end);
end
%%----------------------------------------------------------------------
%% LOOP OVER EACH DATA FILE --------------------------------------------
%%----------------------------------------------------------------------
for i=1:nDataTypes
  %%FIND RUN NAME FOR THE CURRENT DATA FILE
  remains = strtok(dataFiles{i}, '.');
  while length(remains)>0
    [runName, remains] = strtok(remains, '/');
  end
  %%READ IN DATA; SORT BY CLUSTERING PARTITION
  working      = importdata(dataFiles{i}, ',',1);
  data         = working.data;
  data         = data(sortIndex, :);
  if isfield(working, 'featureNames')
    currentNames = working.featureNames;
  else
    currentNames = working.textdata(1,2:end);
  end
  %%BASIC SORT OF FEATUREs
  [dum, index] = sort(mean(data));
  data         = data(:,index);
  currentNames = currentNames(index);
  %%OPTION TO SORT FEATURES USING BIOMARKER PROBS
  if length(featureProbs)>0
    nFeatures    = length(currentNames);
    currentIndex = zeros(1, nFeatures);
    for j=1:nFeatures
      currentIndex(j) = find(strcmp(featureNames, currentNames{j}));
    end
    currentProbs            = featureProbs(currentIndex);
    [currentProbs newIndex] = sort(currentProbs, 'descend');
    data                    = data(:, newIndex);
    currentNames            = currentNames(newIndex);
    %%REMOVE FEATURES THAT ARE LARGELY 'OFF'
    keep         = currentProbs>0.1;
    currentProbs = currentProbs(keep);
    data         = data(:, keep);
    currentNames = currentNames(keep);  
    %%FOR PLOTTING PURPOSES, CLIP OUTLIER VALUES
    lower = quantile(data(:), 0.05);
    upper = quantile(data(:), 0.95);
    data  = max(data, lower);
    data  = min(data, upper);
  end  
  %%PLOT DATA FIGURE
  subplot(dimensionSize, dimensionSize, i)
  imagesc(data'), colorbar
  ylabel('Features')
  xlabel('Items')
  title(runName)
  hold all
  %%LABEL THE CLUSTERING PARTITION
  index                  = find(mod(clusterIDs, 2));
  partitionLabels        = repmat({'.'}, 1, nDataItems);
  partitionLabels(index) = {'l'};
  set(gca, 'XTick',      1:nDataItems)
  set(gca, 'XTickLabel', partitionLabels)
  xlabel('(partition)')
  %%OPTION TO LABEL THE FEATURES
  if length(featureProbs)>0
    featureLabels = cell(1, nFeatures);
    for j=1:10
      threshold    = 0.1*j;
      currentIndex = find(currentProbs>=threshold);
      featureLabels(currentIndex) = strcat(featureLabels(currentIndex), '-');
    end
    set(gca, 'YTick',      1:nFeatures)
    set(gca, 'YTickLabel', featureLabels)
  end  
end
%%----------------------------------------------------------------------
%% SAVE TO FILE --------------------------------------------------------
%%----------------------------------------------------------------------
saveas(gcf, outputFile, 'pdf')
subplot(1,1,1)
hold off
%%*****************************************************************************
%%*** END OF THE FUNCTION *****************************************************
%%*****************************************************************************
%%----------------------------------------------------------------------
%% ----------------------------------------
%%----------------------------------------------------------------------



  %%SORT FEATURES INTO A USEFUL ORDER
%  nFeatures    = size(data, 2);
%  nClusters    = 2*ceil(log(nFeatures));
%  T            = clusterdata(data','maxclust',nClusters); 
%  [dum, index] = sort(T);
%  data         = data(:,index);
