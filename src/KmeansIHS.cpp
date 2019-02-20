#include "RcppArmadillo.h"
using namespace Rcpp;

#include <cmath>
#include <stdlib.h>
#define SQ(x) ((x) * (x))

// [[Rcpp::depends(RcppArmadillo)]]

// Rcpp::sourceCpp("/Users/Nick/Documents/Dissertation/kmeansIHS/KmeansIHS.cpp")

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]] 
double logjac_IHS(NumericMatrix dat, NumericMatrix lambda, IntegerVector cluster, int lambdaType)
{
  int p = lambda.ncol();
  int n = dat.nrow();
  //int K = lambda.nrow();
  NumericVector jacs(p, 0.0);
  NumericVector temp(n, 0.0);
  double lambdaSquared = 0.0;
  NumericVector colSqu(n, 0.0);
  double valSqu = 0;

  if(lambdaType == 1)
  {
    for(int i = 0; i < p; i++)
    {
      colSqu = SQ(dat(_, i));
      lambdaSquared = SQ(lambda(0,i));
      temp = lambdaSquared * colSqu + 1.0;
      temp = log(temp);
      jacs(i) = sum(temp);
    }
  }
  else if(lambdaType == 2)
  {
    for(int i = 0; i < n; i++)
    {
      for(int m = 0; m < p; m++)
      {
        valSqu = SQ(dat(i, m));
        lambdaSquared = SQ(lambda(cluster(i),m));
        temp = lambdaSquared * valSqu + 1.0;
        temp = log(temp);
        jacs(m) += sum(temp);
      }
    }
  }
  
  return (-.5)*sum(jacs);
}

// [[Rcpp::export]] 
double logjac_IHS1(double x, double lambda)
{
  double temp = 0;
  double x_squared = SQ(x);
  double lambdaSquared = SQ(lambda);
  
  temp = lambdaSquared * x_squared + 1.0;
  temp = log(temp);
  
  return (-.5)*temp;
}

// [[Rcpp::export]]
NumericVector IHS(NumericVector col, double lambda)
{
  if(lambda != 0)
  {
    NumericVector res = SQ(lambda)*SQ(col) + 1;
    res = lambda*col + sqrt(res);
    res = log(res)/lambda;
    return res;
  }
  else
  {
    return col;
  }
}

// [[Rcpp::export]]
double IHS1(double val, double lambda)
{
  if(lambda != 0)
  {
    double res = SQ(lambda)*SQ(val) + 1;
    res = lambda*val + sqrt(res);
    res = log(res)/lambda;
    return res;
  }
  else
  {
    return val;
  }
}

// [[Rcpp::export]] 
NumericMatrix IHSAll(NumericMatrix X, NumericMatrix lambda, int lambdaType, IntegerVector cluster)
{
  NumericMatrix X2(X.nrow(), X.ncol());// = clone(X);
  if(lambdaType == 1)
  {
    for(int i = 0; i < lambda.ncol(); i++)
    {
      X2(_,i) = IHS(X(_,i), lambda(0,i));
    }
  }
  else if(lambdaType == 2)
  {
    for(int i = 0; i < X.nrow(); i++)
    {
      for(int j = 0; j < X.ncol(); j++)
      {
        X2(i,j) = IHS1(X(i,j), lambda(cluster(i),j));
      }
    }
  }
  
  return X2;
}

//[[Rcpp::export]]
NumericVector IHS_WSS(NumericMatrix Xihs, IntegerVector clusters, int K, NumericMatrix means)
{
  int n = Xihs.nrow();
  int p = Xihs.ncol();

  NumericVector res(p, 0.0);

  for(int i = 0; i < n; i++)
  {
    for(int j = 0; j < p; j++)
    {
      res(j) += SQ(Xihs(i,j) - means(clusters(i), j));
    }
  }
  return res;
}

//[[Rcpp::export]]
double IHS_WSS_Update1(NumericVector x, IntegerVector clusters, int K, double lambda)
{
  int n = x.length();
  NumericVector x2(n);
  x2 = IHS(x, lambda);
  NumericVector means(K, 0.0);
  double res = 0.0;
  IntegerVector clustCounts(K, 0);
  
  for(int i = 0; i < n; i++)
  {
      means(clusters(i)) += x2(i);
      clustCounts(clusters(i))++;
  }
  
  for(int i = 0; i < K; i++)
  {
    means(i) /= double(clustCounts(i));
  }

  for(int i = 0; i < n; i++)
  {
     res += SQ(x2(i) - means(clusters(i)));
  }
  
  return res;
}

double IHS_WSS_Update(NumericVector x, IntegerVector clusters, int K, NumericVector lambda)
{
  int n = x.length();
  NumericVector x2(n);
  
  for(int i = 0; i < n; i++)
  {
      x2(i) = IHS1(x(i), lambda(clusters(i)));
  }
  
  NumericVector means(K, 0.0);
  double res = 0.0;
  IntegerVector clustCounts(K, 0);
  
  for(int i = 0; i < n; i++)
  {
    means(clusters(i)) += x2(i);
    clustCounts(clusters(i))++;
  }
  
  for(int i = 0; i < K; i++)
  {
    means(i) /= double(clustCounts(i));
  }
  
  for(int i = 0; i < n; i++)
  {
    res += SQ(x2(i) - means(clusters(i)));
  }
  return res;
}

// [[Rcpp::export]]
NumericMatrix T_means(NumericMatrix Xihs, IntegerVector clusters, int k)
{
  int n = Xihs.nrow();
  int p = Xihs.ncol();
  
  // NumericVector temp(p);
  NumericMatrix means(k, p);
  IntegerVector clusterSize(k, 0);

  //  centersihs = IHSAll(centers, lambda);
  
  for(int ii = 0; ii < n; ii++)
  {
    for(int jj = 0; jj < p; jj++)
    {
      means(clusters(ii), jj) += Xihs(ii, jj);
    }
    clusterSize(clusters(ii))++;
  }
  for(int ii = 0; ii < k; ii++)
  {
    means(ii, _) = means(ii, _) / double(clusterSize(ii));
  }
  return means;
}

// [[Rcpp::export]]
NumericMatrix IHS_Distance_Matrix(NumericMatrix Xihs, IntegerVector clusters, NumericMatrix means, int k, NumericMatrix lambda, NumericMatrix X, int lambdaType)
{
  int n = Xihs.nrow();
  int p = Xihs.ncol();

  NumericMatrix dists(n, k);

  for(int ii = 0; ii < n; ii++)
  {
    for(int jj = 0; jj < k; jj++)
    {
      for(int kk = 0; kk < p; kk++)
      {
        dists(ii, jj) += SQ(Xihs(ii, kk) - means(jj, kk));
        /*if(lambdaType == 2)
        {
          dists(ii,jj) += SQ(IHS1(X(ii,kk), lambda(jj,kk)) - means(jj,kk));
        }
        else if(lambdaType == 1)
        {
          dists(ii,jj) += SQ(IHS1(X(ii,kk), lambda(0,kk)) - means(jj,kk));
        }*/
      }
    }
    dists(ii, _) = pow(dists(ii, _), .5);
    //dists(ii, _) = -.5*log(pow(dists(ii, _), .5));
  }
  
  return dists;
}

// [[Rcpp::export]]
arma::mat chooseLambdaMove(arma::cube scores, int lambdaStepType)
{
  int p = scores.n_cols;
  int n = scores.n_rows;
  int slices = scores.n_slices;
  arma::mat moves = arma::zeros<arma::mat>(n, p);

  if(lambdaStepType == 1)
  {
    arma::uvec tempIndex(3);
 //   arma::vec vecScores = scores(arma::span::all, arma::span::all, arma::span::all);
    arma::vec vecScores = arma::vec(scores.memptr(), scores.size());
    arma::uvec order = arma::sort_index(vecScores);
    
    tempIndex = ind2sub(size(scores), order(order.n_elem - 1));
    if(tempIndex(2) == 0)
    {
      moves(tempIndex(0), tempIndex(1))--;
    }
    else if(tempIndex(2) == 2)
    {
      moves(tempIndex(0), tempIndex(1))++;
    }
  }
  else if(lambdaStepType == 2)
  {
    arma::uvec tempIndex(2);
    arma::uvec order(n*slices);
    arma::vec vecScores(n*slices);
    arma::cube subCube(n, 1, slices);
    for(int i = 0; i < p; i++)
    {
      subCube = scores(arma::span::all, arma::span(i, i), arma::span::all);
      vecScores = arma::vec(subCube.memptr(), subCube.size());
      order = arma::sort_index(vecScores);

      tempIndex = ind2sub(size(subCube), order(order.n_elem - 1));

      if(tempIndex(2) == 0)
      {
        moves(tempIndex(0), i)--;
      }
      else if(tempIndex(2) == 2)
      {
        moves(tempIndex(0), i)++;
      }
    }
  }
  else if(lambdaStepType == 3)
  {
    arma::uvec order(slices);
    arma::vec vecScores(slices);
    
    for(int i = 0; i < p; i++)
    {
      for(int j = 0; j < n; j++)
      {
        vecScores = scores.tube(j,i);
        order = arma::sort_index(vecScores);
        
        if(order(order.n_elem-1) == 0)
          moves(j,i)--;
        else if(order(order.n_elem-1) == 2)
          moves(j,i)++;
      }
    }
  }
  return(moves);
}

// [[Rcpp::export]]
IntegerVector getClosest(NumericMatrix dists)
{
  int n = dists.nrow();
 // int p = dists.ncol();
  IntegerVector closest(n);
  
  for(int i = 0; i < n; i++)
  {
    closest(i) = which_min(dists(i,_));
  }
  
  return closest;
}

int which_equal(NumericVector vect, double val)
{
  int res = -1;
  for(int i = 0; i < vect.size(); i++)
  {
    if(vect(i) == val)
    {
      res = i;
      break;
    }
  }
  return(res);
}

int which_equal(IntegerVector vect, int val)
{
  int res = -1;
  for(int i = 0; i < vect.size(); i++)
  {
    if(vect(i) == val)
    {
      res = i;
      break;
    }
  }
  return(res);
}

arma::cube getStepResults(NumericMatrix X, IntegerVector closest, NumericVector wss, NumericMatrix lambda, NumericVector lambdaSeq, IntegerMatrix lambdaIndex, int lambdaType)
{
  int K = lambda.nrow();
  int p = lambda.ncol();
  int n = X.nrow();
  arma::cube res(lambda.nrow(), lambda.ncol(), 3);
  res.fill(-1*DBL_MAX);
  double WSS_i_placeholder;
  double WSS_temp = 0.0;
  NumericVector temp(n);
  if(lambdaType == 2)
  {
    NumericVector lambdaTemp(K);
    for(int i = 0; i < lambda.nrow(); i++)
    {
      for(int j = 0; j < lambda.ncol(); j++)
      {
        if(lambdaIndex(i,j) > 0)
        {
          WSS_i_placeholder = wss(j);
  
          lambdaTemp = lambda(_,j);
          lambdaTemp(i) = lambdaSeq(lambdaIndex(i,j) - 1);
          
          WSS_temp = IHS_WSS_Update(X(_,j), closest, K, lambdaTemp);
          wss(j) = WSS_temp;
          
          lambda(i,j) = lambdaSeq(lambdaIndex(i,j) - 1);
          res(i,j,0) = -1*(n*p/2.0)*log(sum(wss)) + logjac_IHS(X, lambda, closest, lambdaType);
          lambda(i,j) = lambdaSeq(lambdaIndex(i,j));
          
          wss(j) = WSS_i_placeholder;
        }
        
        res(i,j,1) = -1*(n*p/2.0)*log(sum(wss)) + logjac_IHS(X, lambda, closest, lambdaType);
          
        if(lambdaIndex(i,j) < lambdaSeq.size() - 1)
        {
          WSS_i_placeholder = wss(j);
          
          lambdaTemp = lambda(_,j);
          lambdaTemp(i) = lambdaSeq(lambdaIndex(i,j) + 1);
          
          WSS_temp = IHS_WSS_Update(X(_,j), closest, K, lambdaTemp);
          wss(j) = WSS_temp;
          
          lambda(i,j) = lambdaSeq(lambdaIndex(i,j) + 1);
          res(i,j,2) = -1*(n*p/2.0)*log(sum(wss)) + logjac_IHS(X, lambda, closest, lambdaType);
          lambda(i,j) = lambdaSeq(lambdaIndex(i,j));
          
          wss(j) = WSS_i_placeholder;
        }
        
        res(i,j,0) = res(i,j,0) - res(i,j,1);
        res(i,j,2) = res(i,j,2) - res(i,j,1);
        res(i,j,1) = 0;
      }
    }
  }
  else
  {
    for(int j = 0; j < lambda.ncol(); j++)
    {
      if(lambdaIndex(0,j) > 0)
      {
        WSS_i_placeholder = wss(j);
        
        WSS_temp = IHS_WSS_Update1(X(_,j), closest, K, lambdaSeq(lambdaIndex(0,j) - 1));
        
        wss(j) = WSS_temp;
        
        lambda(0,j) = lambdaSeq(lambdaIndex(0,j) - 1);
        res(0,j,0) = -1*(n*p/2.0)*log(sum(wss)) + logjac_IHS(X, lambda, closest, lambdaType);
        lambda(0,j) = lambdaSeq(lambdaIndex(0,j));
        
        wss(j) = WSS_i_placeholder;
      }
      
      res(0,j,1) = -1*(n*p/2.0)*log(sum(wss)) + logjac_IHS(X, lambda, closest, lambdaType);
      
      if(lambdaIndex(0,j) < lambdaSeq.size() - 1)
      {
        WSS_i_placeholder = wss(j);
        
        WSS_temp = IHS_WSS_Update1(X(_,j), closest, K, lambdaSeq(lambdaIndex(0,j) + 1));
        wss(j) = WSS_temp;
        
        lambda(0,j) = lambdaSeq(lambdaIndex(0,j) + 1);
        res(0,j,2) = -1*(n*p/2.0)*log(sum(wss)) + logjac_IHS(X, lambda, closest, lambdaType);
        lambda(0,j) = lambdaSeq(lambdaIndex(0,j));
        
        wss(j) = WSS_i_placeholder;
      }
      
      res(0,j,0) = res(0,j,0) - res(0,j,1);
      res(0,j,2) = res(0,j,2) - res(0,j,1);
      res(0,j,1) = 0;
    }
  }
  
  return(res);
}


// [[Rcpp::export]]
List KM_Transformed_Internal(NumericMatrix X, int K, NumericMatrix centersR, int lambdaType, NumericVector lambdaSeq, int lambdaStepType, NumericMatrix lambdaR, int maxiter)
{
  int n = X.nrow();
  int p = X.ncol();
  
  NumericMatrix X_trans(n, p);
  
  NumericMatrix dists(n, K);
  IntegerVector closest(n);
  
  NumericMatrix centers = clone(centersR);
  NumericMatrix lambda = clone(lambdaR);
  
  NumericMatrix tmeans(K, p);
  NumericVector WSSbyCol(p);
  NumericMatrix oldLambda(K,p);
  NumericMatrix oldCenters(K, p);
  IntegerVector oldClosest(n);
  
  IntegerMatrix lambdaIndex(K, p);
  for(int i = 0; i < p; i++)
  {
    if(lambdaType == 2)
    {
      for(int j = 0; j < K; j++)
      {
        lambdaIndex(j,i) = which_equal(lambdaSeq, lambda(j,i));
      }
    }
    else
      lambdaIndex(0,i) = which_equal(lambdaSeq, lambda(0,i));
  }
  
  arma::cube moveMat = arma::zeros<arma::cube>(K, p, 3);
  arma::mat lambdaMove = arma::zeros<arma::mat>(K, p);
  
  int newPoint;
  NumericVector clusterCount(K);

  int iterCount = 0;
  bool thingsChangedFlag = TRUE;
  
  double transWSS = 0;
  // In lambdaType == 2 case: Initial clusters assigned as closest in euclidean distance.
  // No choice since need clusters to transform. Would love better solution, but don't see one.
  // Would be better if nearest clusters by euclidean distance were passed from R, since not every transformation has identity special case?
  if(lambdaType == 2)
  {
    NumericMatrix tempZeros(K, p);
    
    X_trans = IHSAll(X, tempZeros, 1, closest);
    tmeans = IHSAll(centers, tempZeros, 1, closest);
    
    for(int ii = 0; ii < n; ii++)
    {
      for(int jj = 0; jj < K; jj++)
      {
        for(int kk = 0; kk < p; kk++)
        {
          dists(ii, jj) += SQ(X_trans(ii, kk) - tmeans(jj, kk));
        }
      }
      dists(ii, _) = pow(dists(ii, _), .5);
    }
    
    closest = getClosest(dists);
  }
  // Gets Distances between IHS transformed X and each IHS transformed centers;
  X_trans = IHSAll(X, lambda, lambdaType, closest);
  tmeans = IHSAll(centers, lambda, lambdaType, closest);
  
  // Get Initial clusters

  for(int ii = 0; ii < n; ii++)
  {
    for(int jj = 0; jj < K; jj++)
    {
      for(int kk = 0; kk < p; kk++)
      {
        dists(ii, jj) += SQ(X_trans(ii, kk) - tmeans(jj, kk));
      }
    }
    dists(ii, _) = pow(dists(ii, _), .5);
  }
  // Passes dists as reference and modified it in function. -- Can't use dists as expected anymore without changing function.
  closest = getClosest(dists);
  
  WSSbyCol = IHS_WSS(X_trans, closest, K, tmeans);
  do
  {
    if(iterCount % 5 == 0) Rcpp::checkUserInterrupt();
    
    iterCount++;
    oldLambda = clone(lambda);
    oldCenters = clone(centers);
    oldClosest = clone(closest);
    moveMat = getStepResults(X, closest, WSSbyCol, lambda, lambdaSeq, lambdaIndex, lambdaType);
    
    lambdaMove = chooseLambdaMove(moveMat, lambdaStepType);

  //  std::cout << "Lambdas: \n" << lambda << "\n";
  //  std::cout << "LambdaIndex: \n" << lambdaIndex << "\n";
  //  std::cout << "MoveMat: \n" << moveMat << "\n";
  //  std::cout << "Lambda Move: \n" << lambdaMove << "\n";
  //  std::cout << "Closest: \n" << closest << "\n\n\n";
   
    if(lambdaType == 2)
    {
      for(int i = 0; i < K; i++)
      {
        for(int j = 0; j < p; j++)
        {
          if(lambdaMove(i,j) == -1)
          {
            lambda(i,j) = lambdaSeq(lambdaIndex(i,j) - 1);
            lambdaIndex(i,j)--;
          }
          else if(lambdaMove(i,j) == 1)
          {
            lambda(i,j) = lambdaSeq(lambdaIndex(i,j) + 1);
            lambdaIndex(i,j)++;
          }
        }
      }
    }
    else
    {
      for(int j = 0; j < p; j++)
      {
        if(lambdaMove(0,j) == -1)
        {
          lambda(0,j) = lambdaSeq(lambdaIndex(0,j) - 1);
          lambdaIndex(0,j)--;
        }
        else if(lambdaMove(0,j) == 1)
        {
          lambda(0,j) = lambdaSeq(lambdaIndex(0,j) + 1);
          lambdaIndex(0,j)++;
        }
      }
    }
    
    clusterCount = NumericVector(K);
    for(int i = 0; i < n; i++)
    {
      clusterCount(closest(i))++;
    }
    
    while(any(clusterCount == 0).is_true())
    {
      for(int i = 0; i < K; i++)
      {
        if(clusterCount(i) == 0)
        {
          newPoint = sample(n, 1)(0)-1;
          clusterCount(i)++;
          clusterCount(closest(newPoint))--;
          closest(newPoint) = i;
          Rcout << "Resampled a cluster because no values in a group.\n";
        }
      }
      break;
    }
    
    X_trans = IHSAll(X, lambda, lambdaType, closest);
    tmeans = T_means(X_trans, closest, K);
    dists = IHS_Distance_Matrix(X_trans, closest, tmeans, K, lambda, X, lambdaType);
    closest = getClosest(dists);
    
    tmeans = T_means(X_trans, closest, K);
    clusterCount = NumericVector(K);
    for(int i = 0; i < n; i++)
    {
      clusterCount(closest(i))++;
    }
    WSSbyCol = IHS_WSS(X_trans, closest, K, tmeans);
    transWSS = sum(WSSbyCol);
        
        
    if(lambdaType == 1)
    {
      thingsChangedFlag = !(is_true(all(oldLambda(0,_) == lambda(0,_))) && is_true(all(oldClosest == closest)));
    }
    else if(lambdaType == 2)
    {
      thingsChangedFlag = !(is_true(all(oldLambda == lambda)) && is_true(all(oldClosest == closest)));
    }
  } while(thingsChangedFlag && iterCount < maxiter);
  
  double objScore = 0;
  objScore = -1*(n*p/2.0)*log(sum(WSSbyCol)) + logjac_IHS(X, lambda, closest, lambdaType);

  centers = NumericMatrix(K,p);
  for(int i = 0; i < n; i++)
  {
    for(int j = 0; j < p; j++)
    {
      centers(closest(i),j) += X(i, j);
    }
  }
  for(int i = 0; i < K; i++)
  {
    centers(i,_) = centers(i,_) / (double) clusterCount(i);
  }

  closest = closest + 1;
  List retList = List::create(Named("objectiveScore") = objScore,
                        Named("transWSS") = transWSS,
                        Named("cluster") = closest,
                        Named("centers") = centers,
                        Named("lambda") = lambda,
                        Named("iter") = iterCount,
                        Named("size") = clusterCount,
                        Named("lambdaType") = lambdaType,
                        Named("stepType") = lambdaStepType);
  
  return retList;
}
