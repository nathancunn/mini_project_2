%%(24/1/12 Rich Savage)
%%Function to estimate biomarker probabilities for a run of MDI
%%
function GenerateClinicalMetrics(dataFiles, outputDir, clinicalFile)
%%----------------------------------------------------------------------
%% CONSTRUCT A RUN NAME ------------------------------------------------
%%----------------------------------------------------------------------
remains  = strtok(dataFiles, '.');
while any(strcmp(remains, '')==0)
  [shortFileNames, remains] = strtok(remains, '/');
end
runName = shortFileNames{1};
for i=2:length(shortFileNames)
  runName = [runName, '_', shortFileNames{i}];
end
outputPath = [outputDir, runName];
%%----------------------------------------------------------------------
%% FIND USEFUL VALUES --------------------------------------------------
%%----------------------------------------------------------------------
clusteringFile = [outputPath, '_clusteringPartition.csv'];
outputFile     = [outputPath, '_survivalCurves.pdf'];
%%----------------------------------------------------------------------
%% READ IN CLUSTERING PARTITION ----------------------------------------
%%----------------------------------------------------------------------
working    = importdata(clusteringFile);
itemNames  = working.rowheaders;
clusterIDs = working.data;
nDataItems = length(itemNames);
nClusters  = length(unique(clusterIDs));
%%----------------------------------------------------------------------
%% READ IN CLINICAL OUTCOME INFORMATION --------------------------------
%%----------------------------------------------------------------------
working       = importdata(clinicalFile);
outcomeNames  = working.textdata(1,2:end);
itemNames_all = working.textdata(2:end,1);
outcomes      = working.data;
%%EXTRACT THE 'died' VARIABLE
index = strcmp(outcomeNames, 'died');
died  = outcomes(:, index);
%%EXTRACT THE TIME-TO-EVENT INFORMATION
index       = strcmp(outcomeNames, 'timeToEvent');
timeToEvent = outcomes(:, index);
%%----------------------------------------------------------------------
%% FIND THE REQUIRED OUTCOMES, TIMES TO EVENT --------------------------
%%----------------------------------------------------------------------
keep = zeros(1, nDataItems);
for i=1:nDataItems
  keep(i) = find(strcmp(itemNames_all, itemNames{i}));
end
died        = died(keep);
timeToEvent = timeToEvent(keep);
%%----------------------------------------------------------------------
%% PLOT KAPLAN-MEIER SURVIVAL CURVE ------------------------------------
%%----------------------------------------------------------------------
figure(1)
legendString = {};
for i=1:nClusters
  newString    = ['Cluster ', num2str(i)];
  legendString = {legendString{:}, newString};
  index        = find(clusterIDs==i);
  [y, x]       = ecdf(timeToEvent(index), 'censoring', died(index)==0, 'function', 'survivor');
  stairs(x, y, 'LineWidth', 1.5)
  hold all
end
title('Kaplan-Meier survival curves')
legend(legendString)
hold off
saveas(gcf, outputFile, 'pdf')

disp('******We need a log-rank p-value here!!')

%%*****************************************************************************
%%*** END OF THE FUNCTION *****************************************************
%%*****************************************************************************
%%----------------------------------------------------------------------
%% ----------------------------------------
%%----------------------------------------------------------------------
