%%(RSS, 6/1/12)
%%Function to summarise the MDI MCMC chains, producing the
%%following outputs:
%% - consensus clustering partition 
%% - fusion matrix for the different data types
%% - data plots, ordered by clustering partition
%% - MCMC metrics
%% - if required, biomarker probabilities
%%
%%NOTE:  We wish to summarise the results for each data type, and
%also some aggregates over all data types.  We need to think
%carefully about what are the most usfeul visualisations
%%
%%
function SummariseMdiResults(dataFiles, outputDir, clinicalFile)
%%----------------------------------------------------------------------
%% FIND USEFUL VALUES --------------------------------------------------
%%----------------------------------------------------------------------
burnInFraction = 0.25;
nDataTypes     = length(dataFiles);
%%----------------------------------------------------------------------
%% CONSTRUCT A RUN NAME ------------------------------------------------
%%----------------------------------------------------------------------
remains = strtok(dataFiles, '.');
while any(strcmp(remains, '')==0)
  [shortFileNames, remains] = strtok(remains, '/');
end
runName = shortFileNames{1};
for i=2:length(shortFileNames)
  runName = [runName, '_', shortFileNames{i}];
end
outputPath = [outputDir, runName];
%%----------------------------------------------------------------------
%% FIND THE MCMC FILES -------------------------------------------------
%%----------------------------------------------------------------------
working    = dir([outputPath, '*mcmcSamples.csv']);
nMcmcFiles = length(working);
mcmcFiles  = cell(1, nMcmcFiles);
for i=1:nMcmcFiles
  mcmcFiles{i} = [outputDir, working(i).name];  
  disp(mcmcFiles{i})
end
%%----------------------------------------------------------------------
%% CHECK TO SEE IF THERE ARE ANY FEATURE SELECTION FILES ---------------
%%----------------------------------------------------------------------
working         = dir([outputPath, '*featureParameters.csv']);
nBiomarkerFiles = length(working);
biomarkerFiles  = cell(1, nBiomarkerFiles);
for i=1:nBiomarkerFiles
  biomarkerFiles{i} = [outputDir, working(i).name];  
end
%%----------------------------------------------------------------------
%% OPTIONAL ADDITIONAL OUTPUTS -----------------------------------------
%%----------------------------------------------------------------------
if nBiomarkerFiles>0 
  [featureProbs, featureNames] = GenerateBiomarkerProbabilities(outputPath, biomarkerFiles, nDataTypes, burnInFraction);
else
  featureProbs = [];
  featureNames = [];
end
%%----------------------------------------------------------------------
%% GENERATE THE REQUIRED OUTPUTS ---------------------------------------
%%----------------------------------------------------------------------
clusterIDs = GenerateClusteringPartition(outputPath, mcmcFiles, nDataTypes, burnInFraction);
GenerateMcmcDiagnostics(outputPath, mcmcFiles, nDataTypes, burnInFraction);
GenerateMcmcHistograms(outputPath, mcmcFiles, nDataTypes, burnInFraction);
GenerateFusionMatrix(outputPath, mcmcFiles, nDataTypes, burnInFraction);
GenerateDataPlots(outputPath, dataFiles, clusterIDs, featureProbs, featureNames);
%%----------------------------------------------------------------------
%% ALSO GENERATE PARTITIONS FOR EACH DATA TYPE -------------------------
%%----------------------------------------------------------------------
for i=1:nDataTypes
  currentPath = [outputDir, shortFileNames{i}];
  currentIDs  = GenerateClusteringPartition(currentPath, mcmcFiles, nDataTypes, burnInFraction, i);
  GenerateDataPlots(currentPath, dataFiles(i), currentIDs, featureProbs, featureNames);
end
end
%%*****************************************************************************
%%*** END OF THE FUNCTION *****************************************************
%%*****************************************************************************
%%----------------------------------------------------------------------
%% ----------------------------------------
%%----------------------------------------------------------------------


