function [output1 output2] = BagOfWords(input, mode, geneIndex, dataForCurrentGene)
switch mode
    case 'init'
        data              = input.data;
        nGenesOverall     = input.nGenes;
        nFeatures         = input.nFeatures;
        maxNumberOfComponents = input.maxNumberOfComponents;
        sparseMatrix      = zeros(nGenesOverall,nFeatures);
        sparseVector      = false(1,nGenesOverall);
        
        % Define the cluster structure
        clusterStruct(1,maxNumberOfComponents) = struct(...
            'nFeatures', [], ...
            'nGenesOverall', [], ...
            'beta', [], ...
            'sumBeta',[],...
            'nGenes', [], ...
            'logMarginalLikelihood', [],...
            'dataCounts', [], ...
            'logicalGeneIDs', [], ...
            'N', []);
        
        [clusterStruct.nFeatures       ] = deal(nFeatures);
        [clusterStruct.nGenesOverall   ] = deal(nGenesOverall);

        hyperParameters = 0.5 + zeros(1,nFeatures);
        beta         = hyperParameters;
        [clusterStruct.beta]    = deal(beta);
        [clusterStruct.sumBeta] = deal(sum(beta));

        [clusterStruct.logMarginalLikelihood] = deal(0);
       
        [clusterStruct.nGenes] = deal(0);
 
        [clusterStruct.logicalGeneIDs] = deal(sparseVector);
        
        % Initialise Context 1 clusters:
% % %         nStartingClusters = min(ceil(log(nGenesOverall)), maxNumberOfComponents);
% % %         clusterIDs        = random('unid', nStartingClusters, 1, nGenesOverall); %row vector
        D = pdist(data, 'hamming');
        Z = linkage(D, 'average');
        T = cluster(Z, ceil(log(nGenesOverall)));
        clusterIDs        = T'; %row vector
               
        for i = 1:maxNumberOfComponents
            logicalIndices                   = clusterIDs == i;
            indices                          = find(logicalIndices);
            nGenesInCluster                  = length(indices);
            dataInCluster                    = sparseMatrix;
            dataInCluster(indices,:)         = data(logicalIndices,:);
            currentCluster                   = clusterStruct(i);
            currentCluster.nGenes            = nGenesInCluster;
            currentCluster.dataCounts        = sum(dataInCluster);
            currentCluster.logicalGeneIDs    = logicalIndices;
            currentCluster.N                 = nFeatures*nGenesInCluster;
            if (nGenesInCluster > 0)
                currentCluster = BagOfWords(currentCluster, 'marginal');
            end
            clusterStruct(i) = currentCluster;
        end
        

        
        output1 = clusterStruct;
        output2 = clusterIDs;
        
    case 'marginal'
        dataCounts     = input.dataCounts;
        beta      = input.beta;
        sumBeta   = input.sumBeta;
        nTotalCounts = sum(dataCounts);
        
        logMarginalLikelihood = sum(gammaln(dataCounts+beta));
        logMarginalLikelihood = logMarginalLikelihood + gammaln(sumBeta);
        logMarginalLikelihood = logMarginalLikelihood - sum(gammaln(beta));
        logMarginalLikelihood = logMarginalLikelihood - gammaln(sumBeta + nTotalCounts);
        input.logMarginalLikelihood = logMarginalLikelihood;
        output1 = input;
    
    case 'initialiseAuxiliary'
        output1 = input;
        nGenesInCluster    = 1;
        output1.nGenes     = nGenesInCluster;
        output1.N          = input.nFeatures;
        output1.dataCounts = dataForCurrentGene;
        output1.logicalGeneIDs(geneIndex) = true;
        output1 = BagOfWords(output1, 'marginal');
    case 'removeGene'
        output1 = input;

        %%%%%
        nGenesInCluster     = output1.nGenes - 1;
        dataCounts          = output1.dataCounts - dataForCurrentGene;
        logicalGeneIDs = output1.logicalGeneIDs;
        logicalGeneIDs(geneIndex) = false;
        output1.logicalGeneIDs = logicalGeneIDs;
        
        
        output1.nGenes = nGenesInCluster;
        output1.dataCounts        = dataCounts;
        output1.N                 = nGenesInCluster*output1.nFeatures;
        if (nGenesInCluster == 0)
            output1.logMarginalLikelihood = 0;
        else
            output1 = BagOfWords(output1, 'marginal');
        end
    case 'addGene'
        output1 = input;
        initialNGenes = output1.nGenes;
        if(initialNGenes > 0)
            dataCounts          = output1.dataCounts + dataForCurrentGene;
        else
            dataCounts = dataForCurrentGene;
        end
        %%%%%
        nGenesInCluster     = initialNGenes + 1;
        logicalGeneIDs = output1.logicalGeneIDs;
        logicalGeneIDs(geneIndex) = true;
        output1.logicalGeneIDs = logicalGeneIDs;
        
        
        output1.nGenes            = nGenesInCluster;
        output1.dataCounts        = dataCounts;
        output1.N                 = nGenesInCluster*output1.nFeatures;
        output1 = BagOfWords(output1, 'marginal');
        
end

end







