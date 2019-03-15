#' @name tikmeans_KSelect
#' @aliases tikmeans_KSelect
#' @title TiK-means: Transformation-infused K-means K selection algorithm
#' @description Iterates from K_min to K_max and calculates the jump statistic over various values of $\\eta$ (exps). Prints plots, returns matrix of jump statistics and TiK-means results.
#' @param dat Dataset to cluster
#' @param K_min Lower bound for clusters to check. Default 1.
#' @param K_max Upper bound for clusters to check.
#' @param exps Exponents used in jump statistic calculations.
#' @param lambdaSeq Lambda search grid. Passed to `tikmeans` function
#' @param nstart Number of initial points for tikmeans
#' @param maxiter Max number of iterations in tikmeans
#' @param lambdaType Type of lambda values. Must be 1 or 2. 1 is a single lambda per dimension. 2 is a lambda per dimension*cluster
#' @param lambdaStepType Type of step made on lambda
#' @param verbose Whether to print details of initial points
#' @param numCores Number of cores used in parallel computation.
#' @return List of 3 items
#' \itemize{
#' \item K_tikMeans K selected by jump statistic
#' \item resMat Results of distorted objective function. Each row is calculated with a different exponent from exps.
#' \item km_list TiK-means results for K = K_min to K = K_max
#' }
#' @importFrom graphics lines par plot
#' @export 
tikmeans_KSelect = function(dat, K_min = 1, K_max, exps = seq(0, 100, length.out = 1000), lambdaSeq = seq(0, 10, length.out = 51), lambdaStart = NULL, nstart = 100, maxiter = 10000, lambdaType = 1, lambdaStepType = 1, verbose = FALSE, numCores = detectCores() %/% 2)
{
  km_list <- vector("list", K_max - K_min + 1) 
  WSS_list = rep(NA, K_max - K_min + 1)
  obj_list = rep(NA, K_max - K_min + 1)
  logjac_list = rep(NA, K_max - K_min + 1)
  for(i in K_min:K_max)
  {
    km_list[[i - K_min + 1]] = tikmeans(dat, i, nstart = nstart, lambda = lambdaStart, lambdaSeq = lambdaSeq, lambdaType = lambdaType, lambdaStepType = lambdaStepType, verbose = verbose, numCores = numCores)
    if(lambdaType == 1)
      print(km_list[[i - K_min + 1]]$lambda[1,])
    else
      print(km_list[[i - K_min + 1]]$lambda)
    WSS_list[i - K_min + 1] = km_list[[i - K_min + 1]]$transWSS
    obj_list[i - K_min + 1] = km_list[[i - K_min + 1]]$objectiveScore
    logjac_list[i - K_min + 1] = obj_list[i - K_min + 1] + (prod(dim(dat)))/2*log(WSS_list[i - K_min + 1])
  }
  distortion = WSS_list/prod(dim(dat))
  objtion = -1*obj_list/prod(dim(dat))

  resMat = matrix(NA, ncol = K_max - K_min + 1, nrow = length(exps))
  for(i in 1:length(exps))
  {
    jumpTest = (objtion)^(-exps[i])
    resMat[i,] = c(jumpTest[1], jumpTest[-1] - jumpTest[-length(jumpTest)])
  }
  
  kselVect = (K_min:K_max)[apply(resMat, 1, which.max)]
  plot(exps, kselVect, main = "Jump Statistic By Exponent", xlab = "-1*Exponent", ylab = "K")
  
  if(max(kselVect) != K_max)
  {
    warning("Max value of K never reached. Expand exps. If computation was long use existing km_list returned and calculate value on own with new exps.")
    return(list("K_tikMeans" = NA, "resMat" = NULL, "km_list" = km_list))
  }
  
  topVal = min(which(kselVect == max(kselVect)))
  bottomVal = which(kselVect == min(kselVect))
  bottomVal = max(bottomVal[bottomVal < topVal])
  
  K_tikMeans = -1
  if(topVal > bottomVal + 1)
  {
    kselVect_bw = kselVect[(bottomVal+1):(topVal-1)]
    tab_Ksel = table(kselVect_bw)
    K_tikMeans = min(names(tab_Ksel)[which.max(tab_Ksel)])
  }
  else
  {
    warning("No value chosen between K_min and K_max. Returning K_max.")
    K_tikMeans = K_max
  }

  return(list("K_tikMeans" = K_tikMeans, "resMat" = resMat, "km_list" = km_list))
}

#' @name IHS_BackTransform
#' @aliases IHS_BackTransform
#' @title Back transform dataset with IHS
#' @description Converts a dataset to its transformed counterpart using parameters returned from tikmeans.
#' @param dat Dataset to be back transformed. Must be numeric.
#' @param lambda Transformation parameters.
#' @param lambdaType Dimensions of lambda matrix. 1 means p-dimensional lambda vector. 2 means (kxp)-dimensional lambda matrix. 1 is default.
#' @return Matrix containing observations on the transformed scale.
#' \itemize{
#' \item objectiveScore Value of objective function.
#' \item transWSS WSS on transformed scale
#' \item cluster Cluster assignmets
#' \item centers Cluster centers. On pre-transformation scale.
#' \item lambda Estimated lambda values for IHS transformation.
#' \item iter Number of iterations made by TiK-means
#' \item size Cluster sizes
#' \item lambdaType lambdaType specified in function call
#' \item stepType lambdaStepType specified in function call
#' }
#' @export
IHS_BackTransform = function(dat, lambda, lambdaType = 1, cluster = 1)
{
  dat = as.matrix(dat)
  if(lambdaType == 1)
  {
    for(i in 1:dim(dat)[2])
    {
      dat[,i] = IHS(dat[,i], lambda[1,i])
    }
  }
  else if(lambdaType == 2)
  {
    for(i in 1:dim(dat)[2])
    {
      for(j in 1:dim(dat)[1])
      {
        dat[j,i] = IHS(dat[j,i], lambda[cluster[j],i])
      }
    }
  }
  return(dat)
  
}

# Internal function for determining initial values for lambda. For lambdaType 1, completely random. For lambdaType 2, random across dimensions, but the same for each cluster.
lambdaStarting = function(lambdaSeq, lambdaType, p, K)
{
  retVal = matrix(NA, K, p)
  if(lambdaType == 1)
  {
    retVal[1,] = sample(lambdaSeq, p, replace = TRUE)
  }
  else if(lambdaType == 2)
  {
    retVal = matrix(rep(sample(lambdaSeq, p, replace = TRUE), each = K), nrow = K)
  }
  return(retVal)
}

#' @name tikmeans
#' @aliases tikmeans
#' @title TiK-means: Transformation-infused K-means
#' @description Performs K-means clustering and estimates Inverse Hyperbolic Sine transformation for the provided data and K value. Default values are provided, but setting them based on data is common.
#' @param dat Dataset to be clustered. Must be numeric.
#' @param K Number of clusters to estimate.
#' @param centers Initial cluster means. Must have K rows.
#' @param lambdaType Dimensions of lambda matrix. 1 means p-dimensional lambda vector. 2 means (kxp)-dimensional lambda matrix. 1 is default.
#' @param lambdaStepType Type of step made on lambda. 1 is default and suggested.
#' @param lambdaSeq Lambda search grid. Passed to `tikmeans` function
#' @param lambda Initial values of the lambda matrix.
#' @param nstart Number of initial points for tikmeans.
#' @param maxiter Max number of iterations in tikmeans.
#' @param numCores Number of cores used in parallel computation.
#' @param verbose Whether to print details of initial points
#' @return List containing results from TiK-means.
#' \itemize{
#' \item objectiveScore Value of objective function.
#' \item transWSS WSS on transformed scale
#' \item cluster Cluster assignmets
#' \item centers Cluster centers. On pre-transformation scale.
#' \item lambda Estimated lambda values for IHS transformation.
#' \item iter Number of iterations made by TiK-means
#' \item size Cluster sizes
#' \item lambdaType lambdaType specified in function call
#' \item stepType lambdaStepType specified in function call
#' }
#' @importFrom parallel detectCores
#' @import foreach
#' @import doParallel
#' @useDynLib TiKmeans, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @export
tikmeans = function(dat, K, centers = NULL, lambdaType = 1, lambdaSeq = seq(0,10,length.out= 51), lambdaStepType = 1, lambda = NULL, nstart = 100, maxiter = 1000, numCores = detectCores() %/% 2, verbose = FALSE)
{
  if(lambdaType < 1 | lambdaType > 2)
  {
    stop("Lambda type specified incorrectly. Use 1 for p dimensional lambda, and 2 for p*k dimensional.")
  }
  
  if(lambdaStepType > 3 | lambdaType < 1)
  {
    stop("Lambda step type specified incorrectly. Use 1 for a single update per iteration, 2 for an update per p per iteration, and 3 for an update per k*p per iteration.")
  }
  
  if(any(lambdaSeq < 1e-10 & lambdaSeq > 0))
  {
    cat("Making arbitrarily small lambdaSeq values equal to 0 and removing duplicates.")
    lambdaSeq[lambdaSeq<1e-10] = 0
  }
  lambdaSeq = sort(unique(lambdaSeq))
  
  dat = as.matrix(dat)
  
 # require(doParallel)
  registerDoParallel(numCores)
  cat("Using", numCores, "nodes.\n")

  if(!is.null(centers) & !is.null(lambda))
  {
    nstart = 1
    res = KM_Transformed_Internal(dat, K, centers, lambdaType, lambdaSeq, lambdaStepType, matrix(lambda[c(1:K),], nrow = K), maxiter)
    return(res)
  }
  else if(!is.null(centers))
  {
    res = vector("list", nstart)
    res = foreach(i = 1:nstart) %dopar% KM_Transformed_Internal(dat, K, centers, lambdaType, lambdaSeq, lambdaStepType, lambdaStarting(lambdaSeq, lambdaType, dim(dat)[2], K), maxiter)
  }
  else if(!is.null(lambda))
  {
    res = vector("list", nstart)
    res = foreach(i = 1:nstart) %dopar% KM_Transformed_Internal(dat, K, matrix(dat[sample.int(dim(dat)[1],size = K),], byrow = F, ncol = dim(dat)[2]), lambdaType, lambdaSeq, lambdaStepType, matrix(lambda[c(1:K),], nrow = K), maxiter)
  }
  else
  {
    res = vector("list", nstart)
    res = foreach(i = 1:nstart) %dopar% KM_Transformed_Internal(dat, K, matrix(dat[sample.int(dim(dat)[1],size = K),], byrow = F, ncol = dim(dat)[2]), lambdaType, lambdaSeq, lambdaStepType, lambdaStarting(lambdaSeq, lambdaType, dim(dat)[2], K), maxiter)
  }
  
  if(verbose)
  {
    cat("Printing Summaries:\n")
    cat("\tObjective Functions:\n")
    print(summary(sapply(res, "[[", 1)))
    cat("\tLambdas:\n")
    print(summary(t(sapply(res, "[[", 5))))
    cat("\tIterations:\n")
    print(summary(sapply(res, "[[", 6)))
    cat("\tProportion of random starts that hit iteration limit:\n", mean(unlist(sapply(res, "[[", 6)) >= maxiter))
    cat("\n\n")
    print(hist(sapply(res, "[[", 1), main = "Objective Functions"))
  }
  return(res[[which.max(sapply(res, "[[", 1))]])
}
