%%(6/1/12 Rich Savage)
%%Function to generate MCMC diagnostics for a run of MDI.
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
function GenerateMcmcDiagnostics(outputPath, mcmcFiles, nDataTypes, burnInFraction)
  %%----------------------------------------------------------------------
  %% FIND USEFUL VALUES --------------------------------------------------
  %%----------------------------------------------------------------------
  nParams    = nDataTypes + nDataTypes*(nDataTypes-1)/2;  %%the subset of parameters to display
  nMcmcFiles = length(mcmcFiles);
  outputFile = [outputPath, '_mcmcDiagnostics.pdf'];
  %%----------------------------------------------------------------------
  %% CONSTRUCT PARAMETER NAMES -------------------------------------------
  %%----------------------------------------------------------------------
  paramNames = {};
  %%MASS PARAMETERS
  for i=1:nDataTypes
    newName    = ['massParam_' num2str(i)];
    paramNames = {paramNames{:}, newName};
  end
  %%PHI PARAMETERS
  for i=1:(nDataTypes-1)
    for j=(i+1):nDataTypes
    newName    = ['Phi_' num2str(i) '_' num2str(j)];
    paramNames = {paramNames{:}, newName};
    end
  end
  %%----------------------------------------------------------------------
  %% PLOT PARAMETER TIME SERIES ------------------------------------------
  %%----------------------------------------------------------------------
  %%for ease, assume all files have the same column headers
  figure(1)
  index         = 1:nParams;
  dimensionSize = ceil(sqrt(nParams));
  subplot(dimensionSize, dimensionSize, 1)
  for i=1:nMcmcFiles
    data     = ReadInMcmcFile(mcmcFiles{i});
    data     = data(:, index);
    nSamples = size(data,1);
    keep     = ceil(burnInFraction*nSamples):nSamples;
    data     = data(keep, :);
    for j=1:nParams
      %%TIME-SERIES PLOT
      subplot(dimensionSize, dimensionSize, j)
      plot(data(:, j))
      hold all
    end
  end
  %%----------------------------------------------------------------------
  %% TAG THE PLOTS WITH PARAMETER NAMES ----------------------------------
  %%----------------------------------------------------------------------
  for i=1:nParams
      subplot(dimensionSize, dimensionSize, i)
      title(paramNames{i}, 'FontSize', 6, 'Interpreter', 'none')
      hold off
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
