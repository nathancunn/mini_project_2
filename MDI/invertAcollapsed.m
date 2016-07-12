function [collapsedInvertedA logDetA] = invertAcollapsed(collapsedCovarianceMatrix, nGenes, sigmasq)

b = collapsedCovarianceMatrix(1);
a = b + sigmasq;
n = nGenes;

d = (1 - (a/sigmasq))/((n*b) + sigmasq);

collapsedInvertedA = d;

newsigmasq = sigmasq/b;

% We can uncomment the below if we wish to calculate the determinant rather
% than the log-determinant
%detA = (newsigmasq^n + (n*newsigmasq^(n-1)))*b^n;

logDetA = log(newsigmasq + n) + (n-1)*log(newsigmasq) + n*log(b);

end