%%(RSS, 17/4/13)
%%Script to perform an example MDI run
%%
%%----------------------------------------------------------------------
%% SET UP --------------------------------------------------------------
%%----------------------------------------------------------------------
addpath(genpath('~/MatlabCode/MDI/'));
addpath(genpath('~/MatlabCode/Library/'));
%%----------------------------------------------------------------------
%% SET SOME MDI CONTROL PARAMETERS -------------------------------------
%%----------------------------------------------------------------------
uniqueID          = 1;
nComponents       = 200;%%this can actually be smaller than the
                        %number of data items...
nMcmcSamples      = 1e2;%%this should ideally be 1e4 or more
hyperSamplingFreq = 1;
samplingFreq      = 10;
initialise        = true;
outputPath        = '~/ScienceProjects/MDI/Results/';
%%----------------------------------------------------------------------
%% TIME-SERIES SYNTHETIC DATA ------------------------------------------
%%----------------------------------------------------------------------
%dataFiles = {'Data/MDItestdata1.csv',    'Data/MDItestdata2.csv',    'Data/MDItestdata3.csv'};
%dataTypes = {'TimeCourseUnequalSpacing', 'TimeCourseUnequalSpacing', 'TimeCourseUnequalSpacing'};
%%RUN THE ANALYSIS
%MDI(outputPath, nComponents, nMcmcSamples, dataFiles, dataTypes, ...
%    hyperSamplingFreq, samplingFreq, uniqueID, initialise);
%SummariseMdiResults(dataFiles, outputPath);
%%----------------------------------------------------------------------
%% COUNTS SYNTHETIC DATA (BAG-OF-WORDS MODEL) --------------------------
%%----------------------------------------------------------------------
%dataFiles       = {'Data/PoissonTestData1.csv', 'Data/PoissonTestData2.csv'};
%dataTypes       = {'BagOfWords',                'BagOfWords'};
%featureSwitches = [0, 0];
%%RUN THE ANALYSIS
%MDI(outputPath, nComponents, nMcmcSamples, dataFiles, dataTypes, ...
%    hyperSamplingFreq, samplingFreq, uniqueID, initialise, featureSwitches);
%SummariseMdiResults(dataFiles, outputPath);
%%----------------------------------------------------------------------
%% COUNTS SYNTHETIC DATA (NEGATIVE BINOMIAL MODEL) ---------------------
%%----------------------------------------------------------------------
dataFiles       = {'Data/PoissonTestData1.csv', 'Data/PoissonTestData2.csv'};
dataTypes       = {'Poisson',                   'Poisson'};
%dataTypes       = {'NegativeBinomial',          'NegativeBinomial'};
featureSwitches = [1, 1];
%%RUN THE ANALYSIS
MDI(outputPath, nComponents, nMcmcSamples, dataFiles, dataTypes, ...
    hyperSamplingFreq, samplingFreq, uniqueID, initialise, featureSwitches);
SummariseMdiResults(dataFiles, outputPath);







%%*****************************************************************************
%%*** END OF THE SCRIPT *******************************************************
%%*****************************************************************************
%%----------------------------------------------------------------------
%% ----------------------------------------
%%----------------------------------------------------------------------
