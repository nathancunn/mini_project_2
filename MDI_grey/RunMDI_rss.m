%%(RSS, 6/1/12)
%%Function to run the MDI code on various cancer data sets
%%
%%THOUGHTS
%%--------
%% - output plots of sorted data
%%
%% - it would be really nice for MDI to work for the nTypes=1 case
%% - Other main change:  Store the log posterior in the output file!!
%%
function RunMDI_rss(datasetName, uniqueIdentifier, initialise)
%%----------------------------------------------------------------------
%% SET UP --------------------------------------------------------------
%%----------------------------------------------------------------------
addpath(genpath('~/MatlabCode/MDI/'));
addpath(genpath('~/MatlabCode/Library/'));
%%----------------------------------------------------------------------
%% SET A DEFAULT INITIALISATION ----------------------------------------
%%----------------------------------------------------------------------
if exist('initialise', 'var')==0
  initialise = true;
end
%%----------------------------------------------------------------------
%% DEFINE USEFUL VALUES ------------------------------------------------
%%----------------------------------------------------------------------
nComponents            = 50;
nMcmcSamples           = 1e6;
samplingFrequency      = 10;
hyperSamplingFrequency = 1;
outputPath             = '~/ScienceProjects/MDI/Results/';
%%----------------------------------------------------------------------
%% FIND THE INPUT DATA FILE NAMES --------------------------------------
%%----------------------------------------------------------------------
switch datasetName
 case 'Test', 
  dataPath  = '~/MatlabCode/MDI/Data/';
  dataFiles = {'MDItestdata1.csv', 'MDItestdata2.csv', 'MDItestdata3.csv', ...
               'MDItestdata4.csv', 'MDItestdata5.csv', 'MDItestdata6.csv'};
  dataTypes = {'TimeCourseUnequalSpacing', 'TimeCourseUnequalSpacing', 'TimeCourseUnequalSpacing', ...
               'TimeCourseUnequalSpacing', 'TimeCourseUnequalSpacing', 'TimeCourseUnequalSpacing'};
  dataFiles = strcat(dataPath, dataFiles);
  featureSelectionSwitches = zeros(1,6);

 case 'GaussianTest', 
  dataPath  = '~/MatlabCode/MDI/Data/';
  dataFiles = {'MDItestdata1.csv', 'MDItestdata2.csv', 'MDItestdata3.csv', ...
               'MDItestdata4.csv', 'MDItestdata5.csv', 'MDItestdata6.csv'};
  dataTypes = {'Gaussian', 'Gaussian', 'Gaussian', ...
               'Gaussian', 'Gaussian', 'Gaussian'};
  dataFiles = strcat(dataPath, dataFiles);
  featureSelectionSwitches = ones(1,6);

 case 'SmallTest', 
  dataPath  = '~/MatlabCode/MDI/Data/';
  dataFiles = {'MDItestdata1.csv', 'MDItestdata2.csv'};
  dataTypes = {'TimeCourseUnequalSpacing', 'TimeCourseUnequalSpacing'};
  dataFiles = strcat(dataPath, dataFiles);
  featureSelectionSwitches = zeros(1,2);

 case 'SmallGaussianTest', 
  dataPath  = '~/MatlabCode/MDI/Data/';
  dataFiles = {'MDItestdata1.csv', 'MDItestdata2.csv'};
  dataTypes = {'Gaussian', 'Gaussian'};
  dataFiles = strcat(dataPath, dataFiles);
  featureSelectionSwitches = zeros(1,2);

 case 'GaussianTestData', 
  dataPath  = '~/MatlabCode/MDI/Data/';
  dataFiles = {'GaussianTestData1.csv', 'GaussianTestData2.csv'};
  dataTypes = {'Gaussian', 'Gaussian'};
  dataFiles = strcat(dataPath, dataFiles);
  featureSelectionSwitches = ones(1,2);

 case 'ProstateCancer', 
  dataPath     = '~/Data/CambridgeCancerDatasets/ProstateCancerData/';
  dataFiles    = {'ConcommitantFeatureGE.csv', 'ConcommitantFeatureCN.csv'};
  dataTypes    = {'Multinomial', 'Multinomial'};
  dataFiles    = strcat(dataPath, dataFiles);
  featureSelectionSwitches = ones(1,2);
  
 case 'Metabric', 
  dataPath     = '~/Data/CambridgeCancerDatasets/Metabric/MDI/';
  dataFiles    = {'metabricExpressionData.csv', 'metabricCopyNumberData.csv'};
  dataTypes    = {'Gaussian', 'Gaussian'};
  dataFiles    = strcat(dataPath, dataFiles);
  featureSelectionSwitches = ones(1,2);

 case 'MetabricImage', 
  dataPath     = '~/Data/CambridgeCancerDatasets/Metabric/MDI/';
  dataFiles    = {'metabricExpressionData.csv', 'metabricCopyNumberData.csv', 'metabricImageMetaFeatures.csv'};
  dataTypes    = {'Gaussian', 'Gaussian', 'Gaussian'};
  dataFiles    = strcat(dataPath, dataFiles);
  featureSelectionSwitches = ones(1,3);

 case 'MetabricMoments', 
  dataPath     = '~/Data/CambridgeCancerDatasets/Metabric/MDI/';
  dataFiles    = {'metabricImageMedianData.csv', 'metabricImageSigmaData.csv', 'metabricImageSkewnessData.csv'};
  dataTypes    = {'Gaussian', 'Gaussian', 'Gaussian'};
  dataFiles    = strcat(dataPath, dataFiles);
  featureSelectionSwitches = ones(1,3);

  
  
 case 'RangelBreastCancer', 
  dataPath  = '~/Data/RangelCancerData/';
  dataFiles = {'RangelExpressionData.csv', 'RangelCopyNumberData.csv'};
  dataTypes = {'Gaussian', 'Gaussian'};
  dataFiles = strcat(dataPath, dataFiles);
  featureSelectionSwitches = ones(1,2);

 case 'RangelBreastCancerDiscrete', 
  dataPath  = '~/Data/RangelCancerData/';
  dataFiles = {'RangelExpressionDiscreteData.csv', 'RangelCopyNumberDiscreteData.csv'};
  dataTypes = {'Multinomial', 'Multinomial'};
  dataFiles = strcat(dataPath, dataFiles);
  featureSelectionSwitches = ones(1,2);

 case 'KhanData'
  disp('Under construction!')
  
 case 'Glioblastoma'
  dataPath  = '~/Data/Glioblastoma_TCGA/ProcessedDataFiles/';
  dataFiles = {'GlioblastomaCopyNumber_processed.csv',  'GlioblastomaGeneExpression_processed.csv',...
               'GlioblastomaMethylation_processed.csv', 'GlioblastomaMicroRNA_processed.csv'};
  dataTypes = {'Gaussian', 'Gaussian', 'Multinomial', 'Gaussian'};
  dataFiles = strcat(dataPath, dataFiles);
  featureSelectionSwitches = ones(1,4);
  
 otherwise, disp('Error!  datasetName not recognised')
end
%%----------------------------------------------------------------------
%% TELL THE USER WHICH DATA FILES ARE BEING USED -----------------------
%%----------------------------------------------------------------------
disp(' ')
disp('MULTIPLE DATASET INTEGRATION')
disp('----------------------------')
disp('Running the following data files/types:')
for i=1:length(dataFiles)
  disp(['(', num2str(i), ') ', dataFiles{i}, ' (', dataTypes{i}, ')'])
end
disp(' ')
%%----------------------------------------------------------------------
%% RUN THE MDI ANALYSIS ------------------------------------------------
%%----------------------------------------------------------------------
nCycle      = 1e3;
nIterations = ceil(nMcmcSamples / nCycle);
for i=1:nIterations
  initialise = (i==1) & initialise;
  MDI(outputPath, nComponents, nCycle, dataFiles, dataTypes, ...
      hyperSamplingFrequency, samplingFrequency, uniqueIdentifier, initialise, featureSelectionSwitches);
  SummariseMdiResults(dataFiles, outputPath);
end
end
%%*****************************************************************************
%%*** END OF THE FUNCTION *****************************************************
%%*****************************************************************************
%%----------------------------------------------------------------------
%% ----------------------------------------
%%----------------------------------------------------------------------

%module load matlab
%qsub -t 1-25  -v initialise=1,dataName=MetabricImage      RunMDI.pbs
