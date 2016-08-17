%%(6/1/12 Rich Savage)
%%Function to generate MCMC histograms for a run of MDI.
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
function GenerateMcmcHistograms(outputPath, mcmcFiles, nDataTypes, burnInFraction)
  %%----------------------------------------------------------------------
  %% FIND USEFUL VALUES --------------------------------------------------
  %%----------------------------------------------------------------------
  nParamSubset      = nDataTypes + nDataTypes*(nDataTypes-1)/2;  %%the subset of parameters to display
  nMcmcFiles        = length(mcmcFiles);
  outputFile        = [outputPath, '_mcmcHistograms.pdf'];
  singleDatasetFile = [outputPath, '_singleDatasetHistograms.pdf'];
  %%----------------------------------------------------------------------
  %% PLOT PARAMETER HISTOGRAMS -------------------------------------------
  %%----------------------------------------------------------------------
  %%for ease, assume all files have the same column headers
  keptData      = [];
  index         = 1:nParamSubset;
  dimensionSize = ceil(sqrt(nParamSubset));
  for i=1:nMcmcFiles
    [data, paramNames] = ReadInMcmcFile(mcmcFiles{i});
    data               = data(:, index);
    nSamples           = size(data,1);
    keep               = ceil(burnInFraction*nSamples):nSamples;
    data               = data(keep, :);
    keptData           = [keptData; data];
    %%GENERATE THE PER-DATA-SET HISTOGRAM PLOTS
    figure(2)
    for j=1:nParamSubset
      subplot(dimensionSize, dimensionSize, j)
      hold all
      [yValues, xValues] = hist(data(:,j));
      yValues = 0.9 * yValues / max(yValues);
      plot(xValues, yValues)
      title(paramNames{j}, 'FontSize', 6, 'Interpreter', 'none')
    end
  end
  %%SAVE THE PER DATA SET
  saveas(gcf, singleDatasetFile, 'pdf')
  subplot(1,1,1)
  hold off
  %%TAG THE PLOTS WITH PARAMETER NAMES
  figure(1)
  for i=1:nParamSubset
    subplot(dimensionSize, dimensionSize, i)
    hist(keptData(:, i))
    title(paramNames{i}, 'FontSize', 6, 'Interpreter', 'none')
  end
  %%----------------------------------------------------------------------
  %% SAVE FIGURE TO FILES ------------------------------------------------
  %%----------------------------------------------------------------------
  saveas(gcf, outputFile, 'pdf')
  subplot(1,1,1)
%%*****************************************************************************
%%*** END OF THE FUNCTION *****************************************************
%%*****************************************************************************
%%----------------------------------------------------------------------
%% ----------------------------------------
%%----------------------------------------------------------------------
