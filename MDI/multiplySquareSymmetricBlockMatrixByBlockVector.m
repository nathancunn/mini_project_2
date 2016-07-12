function multipliedBlockMatrices = multiplySquareSymmetricBlockMatrixByBlockVector(compactBlockMatrix1, compactBlockVector, blockSize, nRows)
% CompactBlockMatrix1,2 contain the lower diagonal (including diagonal)
% elements of the block matrix

matrix1 = zeros(nRows, nRows);
i = logical(tril(ones(nRows, nRows)));
matrix1(i) = compactBlockMatrix1;
matrix1 = matrix1.';
matrix1(i) = compactBlockMatrix1;

multipliedBlockMatrices = blockSize*matrix1*compactBlockVector;
end