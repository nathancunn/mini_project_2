%%----------------------------------------------------------------------
%% INVERT BLOCK MATRIX AND CALCULATE DETERMINANT------------------------
%%----------------------------------------------------------------------

function [invertedCollapsedBlockMatrix logBlockMatrixDeterminant] = invertCollapsedBlockMatrix(collapsedCovarianceMatrix, nGenes, nTimes, noise)
[inverseAcollapsed logDetA] = invertAcollapsed(collapsedCovarianceMatrix, nGenes, noise);
if(nTimes == 1)
    invertedCollapsedBlockMatrix = inverseAcollapsed;
    logBlockMatrixDeterminant = logDetA;
else
    b = collapsedCovarianceMatrix(2:nTimes);
    ainverseb = ((1/noise) +nGenes*inverseAcollapsed)*b;
    cainverseb = nGenes*ainverseb*b';
    
    ev = [collapsedCovarianceMatrix((nTimes+1):(end)) - cainverseb(logical(tril(ones(nTimes-1))))];
    [inverseE logDetE] = invertCollapsedBlockMatrix(ev, nGenes, (nTimes-1),noise);
    EinverseC = multiplySquareSymmetricBlockMatrixByBlockVector(inverseE, b, nGenes, nTimes-1)+(b/noise);
    EinverseCAinverse = (nGenes*EinverseC.*inverseAcollapsed) + EinverseC/noise;
    AinverseBEinverseCAinverse = nGenes*ainverseb'*EinverseCAinverse;
    iB1 = inverseAcollapsed+ AinverseBEinverseCAinverse;
    iB3 = -EinverseCAinverse;
    iB4 = inverseE;
    invertedCollapsedBlockMatrix = [iB1; iB3; iB4];
    logBlockMatrixDeterminant = logDetA + logDetE;
end
end


