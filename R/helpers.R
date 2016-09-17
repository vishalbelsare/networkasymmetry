## Convert sparseMatrix to tbl_df.
s_to_df <- function(m) {
  m %>% as("sparseMatrix") %>% summary() %>% tbl_df()
}

## Convert data.frame to a sparseMatrix.
df_to_s <- function(tdf) {
#  with(tdf, sparseMatrix(i=i,j=j,x=x,dims=c(max(i),max(j)))) # dims here.
  sparseMatrix(i=tdf$i,j=tdf$j,x=tdf$x,dims=c(max(tdf$i),max(tdf$j))) # dims here.
}

## Convert either a one-dimensional matrix or a vector to a sparse diagonal matrix.
## Useful alternative to perform elementwise multiplication on large sparse matrices.
to_sdiag <- function(x) { 
  if (class(x)=="dgCMatrix") {
    # If x is a matrix:
    return(.sparseDiagonal(n=length(x[,1]),x=x[,1]))
  } else if (is.vector(x)) {
    # If x is a vector:
    return(.sparseDiagonal(n=length(x),x=x))
  }
}