library(rARPACK)
set.seed(123)
## Some random data
x = matrix(rnorm(1000 * 1000), 1000)
## If retvec == FALSE, we don't calculate eigenvectors
eigs_sym(cov(x), k = 5, which = "LM", opts = list(retvec = FALSE))
#For really large data, the matrix is usually in sparse form. rARPACK
#supports several sparse matrix types defined in the Matrix
#package, and you can even pass an implicit matrix defined by a function to
#eigs(). See ?rARPACK::eigs for details.

library(Matrix)
spmat = as(cov(x), "dgCMatrix")
z <- eigs_sym(spmat, k=1)
str(z)
hist(z$vectors[,1])

## Implicitly define the matrix by a function that calculates A %*% x
## Below represents a diagonal matrix diag(c(1:10))
fmat = function(x, args)
{
  return(x * (1:10))
}
eigs_sym(fmat, 3, n = 10, args = NULL)