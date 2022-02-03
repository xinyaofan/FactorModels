#'borrowed from packaage CopulaModel
#'@description transform data into u-scale
#'@param data input data
#'@param aunif set to be default value -0.5
#'@return out data after transformation
uscore<-function (data, aunif = -0.5)
{
	if (is.vector(data)) {
		nr = length(data)
		us = ((1:nr) + aunif)/(nr + 1 + 2 * aunif)
		jj = rank(data)
		out = us[jj]
	}
	else {
		nc = ncol(data)
		nr = nrow(data)
		out = matrix(0, nr, nc)
		us = ((1:nr) + aunif)/(nr + 1 + 2 * aunif)
		for (j in 1:nc) {
			jj = rank(data[, j])
			tem = us[jj]
			out[, j] = tem
		}
	}
	out
}




isposdef<-function (amat)
{
	tt = try(chol(amat), silent = T)
	ifelse(class(tt) == "matrix", T, F)
}
