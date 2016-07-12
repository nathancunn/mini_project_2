%%(24/1/12 Rich Savage)
%%Function to estimate biomarker probabilities for a run of MDI
%%
function [biomarkerMeans featureNames] = GenerateBiomarkerProbabilities(outputPath, biomarkerFiles, nDataTypes, burnInFraction)
%%----------------------------------------------------------------------
%% FIND USEFUL VALUES --------------------------------------------------
%%----------------------------------------------------------------------
nBiomarkerFiles     = length(biomarkerFiles);
outputFile          = [outputPath, '_biomarkerProbabilities.csv'];
timeseriesFile      = [outputPath, '_nSwitchedOn.pdf'];
[dum, featureNames] = ReadInMcmcFile(biomarkerFiles{1});
nFeatures           = length(featureNames);
probArray           = zeros(nBiomarkerFiles, nFeatures);
nSwitchedOnArray    = zeros(nBiomarkerFiles, nFeatures);
%%----------------------------------------------------------------------
%% ESTIMATE BIOMARKER PROBABILITIES ------------------------------------
%%----------------------------------------------------------------------
%%do this in turn for each input file
%%find mean, std error for each probability
for i=1:nBiomarkerFiles
  data           = ReadInMcmcFile(biomarkerFiles{i});
  nSamples       = size(data,1);
  keep           = ceil(burnInFraction*nSamples):nSamples;
  data           = data(keep, :);
  probArray(i,:) = mean(data);
  newSwitchedOn  = sum(data, 2);
  %%TIME-SERIES PLOT
  plot(newSwitchedOn)
  ylim([0, 1.1*nFeatures])
  title('nSwitchedOn')
  xlabel('MCMC sample')
  ylabel('nFeaturesSwitchedOn')
  hold all
end
biomarkerMeans  = mean(probArray, 1);
if nBiomarkerFiles>1 
  biomarkerSigmas = std(probArray, 1);
else
  biomarkerSigmas = zeros(nFeatures, 1);
end
%%----------------------------------------------------------------------
%% SAVE THE TIME SERIES PLOT TO FILE -----------------------------------
%%----------------------------------------------------------------------
saveas(gcf, timeseriesFile, 'pdf')
hold off
%%----------------------------------------------------------------------
%% SORT THE FEATURES BY BIOMARKER PROBABILITY --------------------------
%%----------------------------------------------------------------------
[biomarkerMeans, sortIndex] = sort(biomarkerMeans', 'descend');
biomarkerSigmas             = biomarkerSigmas(sortIndex);
featureNames                = featureNames(sortIndex);
%%----------------------------------------------------------------------
%% WRITE THE BIOMARKER PROBABILITIES TO FILE ---------------------------
%%----------------------------------------------------------------------
fid = fopen(outputFile, 'wt');
for i=1:nFeatures
  fprintf(fid, [featureNames{i}, ',', num2str(biomarkerMeans(i)), ',', num2str(biomarkerSigmas(i)), '\n']); 
end
fclose(fid);
%%*****************************************************************************
%%*** END OF THE FUNCTION *****************************************************
%%*****************************************************************************
%%----------------------------------------------------------------------
%% ----------------------------------------
%%----------------------------------------------------------------------
