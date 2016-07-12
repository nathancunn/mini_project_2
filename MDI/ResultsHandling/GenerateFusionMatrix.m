%%(6/1/12 Rich Savage)
%%Function to generate the fusion matrix for and MDI run.
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
function GenerateFusionMatrix(outputPath, mcmcFiles, nDataTypes, burnInFraction)
%%----------------------------------------------------------------------
%% FIND USEFUL VALUES --------------------------------------------------
%%----------------------------------------------------------------------
nParams      = nDataTypes + nDataTypes*(nDataTypes-1)/2;  
fusionMatrix = zeros(nDataTypes, nDataTypes);
nMcmcFiles   = length(mcmcFiles);
outputFile   = [outputPath, '_fusionMatrix.pdf'];
%%----------------------------------------------------------------------
%% READ IN THE PHI PARAMETERS ------------------------------------------
%%----------------------------------------------------------------------
%%for ease, assume all files have the same column headers
keptData = [];
index    = (nDataTypes+1):nParams;
for i=1:nMcmcFiles
  data     = ReadInMcmcFile(mcmcFiles{i});
  data     = data(:, index);
  nSamples = size(data,1);
  keep     = ceil(burnInFraction*nSamples):nSamples;
  data     = data(keep, :);
  keptData = [keptData; data];
end
%%----------------------------------------------------------------------
%% CONSTRUCT THE FUSION MATRIX -----------------------------------------
%%----------------------------------------------------------------------
counter = 0;
for i=1:(nDataTypes-1)
  for j=(i+1):nDataTypes
    counter           = counter+1;
    newElement        = mean(data(:,counter));
    fusionMatrix(i,j) = newElement;
    fusionMatrix(j,i) = newElement;
  end
end
%%FILL IN THE DIAGONAL ELEMENTS
%%(Phi is not defined here, so just pick a convenient value for plotting)
maxValue = ceil(1.1 * max(fusionMatrix(:)));
for i=1:nDataTypes
  fusionMatrix(i,i) = maxValue;
end
%%----------------------------------------------------------------------
%% GENERATE, SAVE PLOT -------------------------------------------------
%%----------------------------------------------------------------------
figure(1)
imagesc(fusionMatrix)
colorbar
title('Fusion Matrix')
saveas(gcf, outputFile, 'pdf')
%%*****************************************************************************
%%*** END OF THE FUNCTION *****************************************************
%%*****************************************************************************
%%----------------------------------------------------------------------
%% ----------------------------------------
%%----------------------------------------------------------------------
