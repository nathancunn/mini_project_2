function [output1 output2 output3] = TimeCourse(input, mode, geneIndex, dataForCurrentGene, newLogHypers, newCovarianceMatrixInverses)
switch mode
    case 'init'
        data              = input.data;
        nGenes            = input.nGenes;
        nFeatures         = input.nFeatures;
        sparseMatrix      = zeros(nGenes,nFeatures);
        sparseVector      = false(1,nGenes);
        maxNumberOfComponents = input.maxNumberOfComponents;
        
        featureNames      = input.featureNames;
        featureNames      = cellfun(@str2num,featureNames);
        [X, Y] = meshgrid(featureNames);
        timeDiffs         = (-(X - Y).^2);
        % We assume independent normal priors for the log hyperparameters
        
        % Note: order of hypers: length scale, signal variance,
        % noise variance
        hyperPriorParameters = [0, 1; 0, 1; 0, 1]; % [mean s.d.; ...]
        lowerTriangularLogicalMatrix = logical(tril(ones(nFeatures)));
        % Define the cluster structure
        clusterStruct(1,(maxNumberOfComponents+1)) = struct(...
            'nFeatures', [], ...
            'nGenesOverall', [], ...
            'timeDiffs', [],...
            'logHypers', [], ...
            'logPriorOfLogHypers', [], ...
            'squaredHypers', [], ...
            'hyperPriorParams', [], ...
            'lowerTriangularPartOfCovarianceMatrix', [], ...
            'covarianceMatrixInverses', [], ...
            'nGenes', [], ...
            'logMarginalLikelihood', [],...
            'dataCounts', [], ...
            'squaredDataCounts', [], ...
            'logicalGeneIDs', [], ...
            'lowerTriangularLogicalMatrix', [], ...
            'N', []);
        
        [clusterStruct.nFeatures       ] = deal(nFeatures);
        [clusterStruct.nGenesOverall   ] = deal(nGenes);
        [clusterStruct.hyperPriorParams] = deal(hyperPriorParameters);
        [clusterStruct.timeDiffs] = deal(timeDiffs);
        [clusterStruct.lowerTriangularLogicalMatrix] = deal(lowerTriangularLogicalMatrix);
        
        [clusterStruct.logMarginalLikelihood] = deal(0);
        
        [clusterStruct.nGenes] = deal(0);
        [clusterStruct.logicalGeneIDs] = deal(sparseVector);
        
        % Initialise clusters:
        nStartingClusters = ceil(log(nGenes));
        clusterIDs        = random('unid', nStartingClusters, 1, nGenes); %row vector
        uniqueIDs         = unique(clusterIDs);
        
        for i = 1:maxNumberOfComponents
            clusterStruct(i).covarianceMatrixInverses(1,nGenes) =...
                struct('invertedCovarianceMatrix', [], 'determinant', []);
        end
        
        
        for i = uniqueIDs
            logicalIndices                   = clusterIDs == i;
            indices                          = find(logicalIndices);
            nGenesInCluster                  = length(indices);
            dataInCluster                    = sparseMatrix;
            dataInCluster(indices,:)         = data(logicalIndices,:);
            currentCluster                   = clusterStruct(i);
            currentCluster.logicalGeneIDs    = logicalIndices;
            currentCluster.dataCounts        = sum(dataInCluster,1);
            currentCluster.squaredDataCounts = sum(dataInCluster.^2,1);
            currentCluster.nGenes            = nGenesInCluster;
            currentCluster.N                 = nFeatures*nGenesInCluster;
            
            logHypers = TimeCourse(currentCluster, 'sampleHypers');
            currentCluster.logHypers           = logHypers;
            currentCluster.logPriorOfLogHypers = [log(normpdf(logHypers(1), hyperPriorParameters(1,1), hyperPriorParameters(1,2))),...
                log(normpdf(logHypers(2), hyperPriorParameters(2,1), hyperPriorParameters(2,2))),...
                log(normpdf(logHypers(3), hyperPriorParameters(3,1), hyperPriorParameters(3,2)))];
            currentCluster.squaredHypers = exp(2*currentCluster.logHypers);
            
            
            l2  = currentCluster.squaredHypers(1);
            sf2 = currentCluster.squaredHypers(2);
            
            covarianceMatrix = sf2*exp(timeDiffs/(2*l2));
            lowerTriangularPart = covarianceMatrix(lowerTriangularLogicalMatrix);
            currentCluster.lowerTriangularPartOfCovarianceMatrix = lowerTriangularPart;
            
            currentCluster.covarianceMatrixInverses(1,nGenes) =...
                struct('invertedCovarianceMatrix', [], 'determinant', []);
            
            currentCluster = TimeCourse(currentCluster, 'invert');
            currentCluster = TimeCourse(currentCluster, 'marginal');
            clusterStruct(i) = currentCluster;
        end
        clusterStruct(end) = [];
        
        output1 = clusterStruct;
        output2 = clusterIDs;
    case 'marginal'
        nTimes             = input.nFeatures;
        nGenesInCluster    = input.nGenes;
        se2                = input.squaredHypers(3);
        N                  = input.N;
        dataCounter        = input.dataCounts;
        dataSquaredCounter = input.squaredDataCounts;
        L                  = input.covarianceMatrixInverses(nGenesInCluster).invertedCovarianceMatrix;
        invertedMatrix     = zeros(nTimes, nTimes);
        lowerTriangularLogicalMatrix                 = input.lowerTriangularLogicalMatrix;
        invertedMatrix(lowerTriangularLogicalMatrix) = L;
        invertedMatrix     = invertedMatrix.';
        invertedMatrix(lowerTriangularLogicalMatrix) = L;
        logDetK            = input.covarianceMatrixInverses(nGenesInCluster).determinant;
        yColSum            = dataCounter';
        y2                 = sum(dataSquaredCounter);
        input.logMarginalLikelihood = -0.5*(N*log(2*pi)+logDetK+(yColSum'*invertedMatrix*yColSum + y2/se2));
        output1  = input;
        output2  = [];
    case 'invert'
        nGenesInCluster    = input.nGenes;
        se2                = input.squaredHypers(3);
        nTimes             = input.nFeatures;
        collapsedCovarianceMatrix = input.lowerTriangularPartOfCovarianceMatrix;
        
        
        [invertedK logDetK] =...
            invertCollapsedBlockMatrix(collapsedCovarianceMatrix, nGenesInCluster, nTimes, se2);
        % % % %         %Uncomment the below to check that the inversion is working
        % % % %         %Un-collapse the collapsedCovarianceMatrix
        % % % %         logicalMat = logical(tril(ones(nTimes)));
        % % % %         matrix1 = zeros(nTimes);
        % % % %         matrix1(logicalMat) = collapsedCovarianceMatrix;
        % % % %         matrix1 = matrix1.';
        % % % %         matrix1(logicalMat) = collapsedCovarianceMatrix;
        % % % %
        % % % %         K = matrix1;
        % % % %         fullCovMatrix = blockify(K, nGenesInCluster, false);
        % % % %         invC = inv(fullCovMatrix + se2*eye(size(fullCovMatrix,1)));
        % % % %         invCreduced = invC(2:nGenesInCluster:end, 1:nGenesInCluster:end);
        % % % %
        % % % %         %Un-collapse the matrix
        % % % %         logicalMat = logical(tril(ones(nTimes)));
        % % % %         matrix1 = zeros(nTimes);
        % % % %         matrix1(logicalMat) = invertedK;
        % % % %         matrix1 = matrix1.';
        % % % %         matrix1(logicalMat) = invertedK;
        % % % %         %Should have matrix1 == invCreduced(to many decimal places)
        
        if(~isreal(logDetK))  %We should not ever enter here, but just in case
            disp('Numerical error - covariance matrix may not be positive definite')
        end
        input.covarianceMatrixInverses(nGenesInCluster).invertedCovarianceMatrix...
            = invertedK;
        input.covarianceMatrixInverses(nGenesInCluster).determinant = logDetK;
        output1 = input;
        output2  = [];
    case 'sampleHypers'
        logHypers         = zeros(3,1);
        hyperPriorParameters = input.hyperPriorParams;
        for j  = 1:3
            % Initialise the hypers by sampling from the prior
            hyperPriorParams             = hyperPriorParameters(j,:);
            logHypers(j) ...
                = hyperPriorParams(1) + (hyperPriorParams(2)*randn);
        end
        output1 = logHypers;
        output2 = [];
    case 'removeGene'
        output1           = input;
        nGenesInCluster   = output1.nGenes - 1;
        dataCounts        = output1.dataCounts - dataForCurrentGene;
        squaredDataCounts = output1.squaredDataCounts - dataForCurrentGene.^2;
        logicalGeneIDs = input.logicalGeneIDs;
        logicalGeneIDs(geneIndex) = false;
        output1.logicalGeneIDs = logicalGeneIDs;
        output1.nGenes            = nGenesInCluster;
        output1.dataCounts        = dataCounts;
        output1.squaredDataCounts = squaredDataCounts;
        output1.N                 = nGenesInCluster*output1.nFeatures;
        if (nGenesInCluster == 0)
            output1.logMarginalLikelihood = 0;
        else
            if(isempty(output1.covarianceMatrixInverses(nGenesInCluster).determinant))
                output1 = TimeCourse(output1, 'invert' );
            else
%                 abc = (output1.covarianceMatrixInverses(nGenesInCluster).determinant);
%                 output1 = TimeCourse(output1, 'invert');
%                 if(abs(abc - output1.covarianceMatrixInverses(nGenesInCluster).determinant)>0)
%                     disp('debug!')
%                 end
            end
            output1 = TimeCourse(output1, 'marginal');
        end
    case 'addGene'
        output1 = input;
        initialNGenes = output1.nGenes;
        hyperPriorParameters = output1.hyperPriorParams;
        if(initialNGenes>0)
            dataCounts = output1.dataCounts ...
                + dataForCurrentGene;
            squaredDataCounts = output1.squaredDataCounts ...
                + dataForCurrentGene.^2;
        else
            dataCounts = dataForCurrentGene;
            squaredDataCounts = dataForCurrentGene.^2;
            
            logHypers = TimeCourse(output1, 'sampleHypers');
            output1.logHypers           = logHypers;
            output1.logPriorOfLogHypers = [log(normpdf(logHypers(1), hyperPriorParameters(1,1), hyperPriorParameters(1,2))),...
                log(normpdf(logHypers(2), hyperPriorParameters(2,1), hyperPriorParameters(2,2))),...
                log(normpdf(logHypers(3), hyperPriorParameters(3,1), hyperPriorParameters(3,2)))];
            output1.squaredHypers = exp(2*output1.logHypers);
            l2  = output1.squaredHypers(1);
            sf2 = output1.squaredHypers(2);
            
            covarianceMatrix = sf2*exp(output1.timeDiffs/(2*l2));
            lowerTriangularPart = covarianceMatrix(output1.lowerTriangularLogicalMatrix);
            output1.lowerTriangularPartOfCovarianceMatrix = lowerTriangularPart;
            
            
        end
        
        nGenesInCluster   =  initialNGenes + 1;
        logicalGeneIDs = output1.logicalGeneIDs;
        logicalGeneIDs(geneIndex) = true;
        output1.logicalGeneIDs = logicalGeneIDs;
        
        output1.nGenes            = nGenesInCluster;
        output1.dataCounts        = dataCounts;
        output1.squaredDataCounts = squaredDataCounts;
        output1.N                 = nGenesInCluster*output1.nFeatures;
        if(isempty(output1.covarianceMatrixInverses(nGenesInCluster).determinant))
            output1 = TimeCourse(output1, 'invert');
        else
%             abc = (output1.covarianceMatrixInverses(nGenesInCluster).determinant);
%             output1 = TimeCourse(output1, 'invert');
%             if(abs(abc - output1.covarianceMatrixInverses(nGenesInCluster).determinant)>0)
%                 disp('debug!')
%             end
            
        end
        output1 = TimeCourse(output1, 'marginal');
end

end