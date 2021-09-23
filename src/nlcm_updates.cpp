// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <Rmath.h>
using namespace Rcpp;
using namespace arma;
using namespace R;

// utility functions;

//' xlogx
//'
//' utility function
//'
//' @param x a positive number or zero
//'
//' @useDynLib doubletree
//' @importFrom Rcpp sourceCpp
//' @export
// [[Rcpp::export]]
double xlogx(double x){
  double res=0.0;
  if (x!=0.0){
    res = x*log(x);
  }
  return(res);
}

// [[Rcpp::export]]
double logsumexp(arma::vec logv_arma)
{
  //int n = logv.size();
  //if(n<=1)
  //	cout<<"Warning in logsumexp"<<endl;
  double max = logv_arma.max();
  double answer = 0.0;
  // log(sum(exp(logf)) 	= log(sum(exp(logf - max(logf) + max(logf)))
  //			= max(logf) + log(sum(exp(logf - max(logf)))
  answer = max + log(sum(exp(logv_arma-max)));
  return answer;
}

//' logexpit to avoid numerical underflow
//'
//' @param x a number
//' @useDynLib doubletree
//' @importFrom Rcpp sourceCpp
// [[Rcpp::export]]
double logexpit_cpp(double x)
{
  arma::vec tmp(2);tmp.zeros();
  tmp(1) = -x;
  return(-logsumexp(tmp));
}


// [[Rcpp::export]]
double logsumexp_row(arma::rowvec logv_arma)
{
  //int n = logv.size();
  //if(n<=1)
  //	cout<<"Warning in logsumexp"<<endl;
  double max = logv_arma.max();
  double answer = 0.0;
  // log(sum(exp(logf)) 	= log(sum(exp(logf - max(logf) + max(logf)))
  //			= max(logf) + log(sum(exp(logf - max(logf)))
  answer = max + log(sum(exp(logv_arma-max)));
  return answer;
}

// use the functionality described at http://arma.sourceforge.net/docs.html#each_colrow
// to "Apply a vector operation to each column or row of a matrix "
// also see: https://stackoverflow.com/questions/68351041/how-can-i-replicate-rs-functionality-with-multiplying-a-matrix-by-a-vector-elem
// [[Rcpp::export]]
arma::mat mtv(arma::mat mat, arma::vec v) {
  mat.each_col() %= v;
  return mat;
}

// A <- matrix(c(1,2,3,4,5,6), nrow=2)
//   b <- c(1,3)
//   A * b
//   mtv(A, b)

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

//' Calculate variational moments during the updates
//'
//' Get all moments that need updating when iterating over internal and leaf nodes;
//' for both trees. When updating a node u, all the moments
//' of the descedant nodes will be changed. A recalculation of the moments are necessary
//' when moving on to another node.
//'
//' @param prob1,prob2 variational probabilities; `prob1` is for \code{s*_u} - length `p1`;
//' `prob2` is for `s_cu` - a matrix `pL1` by `p2`; in R, a list of pL1 length - each element being of length `p2`.
//' @param mu_gamma variational Gaussian means (for \code{s*_u=1} component) for J*K
//' logit(class-specific response probabilities); (J,K,p1) array; In R, we used a list of p1 (J,K) matrices
//' @param sigma_gamma variational Gaussian variances (for \code{s*_u=1} component)
//' for J*K logit(class-specific response probabilities); (J,K,p1) array; in R, we used a list o f p1 (J,K) matrices
//' @param mu_alpha variational Gaussian mean vectors (for \code{s_cu=1} component) -
//' this is a pL1 by K-1 by p2 array; in R, we used a list of p2 matrices (each of dimension pL1 by K-1)
//' @param sigma_alpha variational Gaussian variances (for \code{s_cu=1} component)
//' - this is an array of dimension (pL1, K-1, p2); in R, we used a list of p2 matrices,
//' each of dimension pL1 by K-1.
//' @param anc1,anc2 `anc1` is a list of pL1 vectors, each vector has the node ids of the ancestors in tree1;
//' lengths may differ. The ancestors include the node concerned; similarly for `anc2`
//' @param cardanc1,cardanc2 `cardanc1` is a numeric vector of length pL1; integers. The number
//' of ancestors for each leaf node in tree1; similarly for `cardanc2`.
//'
//' @return a List
//'
//' \describe{
//'   return List::create(Named("E_beta")=E_beta,
//'    Named("E_beta_sq")=E_beta_sq,
//'    Named("E_eta")=E_eta,
//'    Named("E_eta_sq")=E_eta_sq);
//'}
//'
//' @useDynLib doubletree
//' @importFrom Rcpp sourceCpp
//' @export
// [[Rcpp::export]]
List get_moments_cpp_doubletree(arma::vec prob1,// p1
                                arma::mat prob2,// pL1 by p2
                                arma::cube mu_gamma,//J by K by p1
                                arma::cube sigma_gamma,//J by K by p1
                                arma::cube mu_alpha, // pL1 by K-1 by p2
                                arma::cube sigma_alpha,//  pL1 by K-1 by p2
                                List anc1, // length pL1; each element is a vector.
                                List anc2, // length pL2; each element is a vector.
                                arma::vec cardanc1, // length pL1; integers
                                arma::vec cardanc2 // length pL2; integers
){
  int pL1 = cardanc1.size();
  int pL2 = cardanc2.size();
  int J = mu_gamma.n_rows;
  int K = mu_gamma.n_cols;
  int p1 = mu_gamma.n_slices;
  int p2 = mu_alpha.n_slices;

  arma::cube E_beta(J,K,pL1);E_beta.zeros();// leaf level
  arma::cube E_beta_sq(J,K,pL1);E_beta_sq.zeros(); // leaf level

  arma::cube  E_eta(pL1,K-1,pL2);E_eta.zeros(); // leaf level
  arma::cube  E_eta_sq(pL1,K-1,pL2);E_eta_sq.zeros();// leaf level

  // follow tree1:
  int n_anc=0;
  int uu=0;
  for (int v=0;v<pL1;v++){
    arma::vec curr_anc = anc1[v];
    n_anc = (int) cardanc1(v);
    for (int u=0;u<n_anc;u++){
      uu =  (int) curr_anc(u)-1;
      E_beta.slice(v)    += prob1(uu)*mu_gamma.slice(uu);
      E_beta_sq.slice(v) += prob1(uu)*(sigma_gamma.slice(uu)+(1.0-prob1(uu))*pow(mu_gamma.slice(uu),2.0)); //not yet.
    }
    E_beta_sq.slice(v)   += pow(E_beta.slice(v),2.0);
  }

  // follow tree2:
  for (int v=0;v<pL2;v++){
    arma::vec curr_anc = anc2[v];
    n_anc = (int) cardanc2(v);
    for (int u=0;u<n_anc;u++){
      uu =  (int) curr_anc(u)-1;
      E_eta.slice(v)        += mtv(mu_alpha.slice(uu),prob2.col(uu));
      E_eta_sq.slice(v)     += mtv(sigma_alpha.slice(uu)+
        mtv(pow(mu_alpha.slice(uu),2.0),(1.0-prob2.col(uu))),
        prob2.col(uu));
    }
    E_eta_sq.slice(v)        += pow(E_eta.slice(v),2.0);
  }
  // return results:
  return List::create(Named("E_beta")=E_beta,
                      Named("E_beta_sq")=E_beta_sq,
                      Named("E_eta")=E_eta,
                      Named("E_eta_sq")=E_eta_sq
  );
}


//' Calculate variational moments during the updates
//'
//' Get all moments that need updating when iterating over internal and leaf nodes;
//' for both trees. When updating a node u, all the moments
//' of the descedant nodes will be changed. A recalculation of the moments are necessary
//' when moving on to another node.
//'
//' @param prob1 variational probabilities; `prob1` is for \code{s*_u} - length `p1`.
//' @param mu_gamma variational Gaussian means (for \code{s*_u=1} component) for J*K
//' logit(class-specific response probabilities); (J,K,p1) array; In R, we used a list of p1 (J,K) matrices
//' @param sigma_gamma variational Gaussian variances (for \code{s*_u=1} component)
//' for J*K logit(class-specific response probabilities); (J,K,p1) array; in R, we used a list o f p1 (J,K) matrices
//' @param E_beta,E_beta_sq J by K by pL1
//' each of dimension pL1 by K-1.
//' @param anc1 `anc1` is a list of pL1 vectors, each vector has the node ids of the ancestors in tree1;
//' lengths may differ. The ancestors include the node concerned; similarly for `anc2`
//' @param leaf_desc an integer vector of leaf descendant ids under consideration
//'
//' @return a List
//'
//' \describe{
//'   return List::create(Named("E_beta")=E_beta,
//'    Named("E_beta_sq")=E_beta_sq);
//'}
//'
//' @useDynLib doubletree
//' @importFrom Rcpp sourceCpp
//' @export
// [[Rcpp::export]]
List get_moments_cpp_eco_gamma_doubletree(arma::vec prob1,// p1
                                          arma::cube mu_gamma,//J by K by p1
                                          arma::cube sigma_gamma,//J by K by p1
                                          arma::cube E_beta,//J by K by pL1
                                          arma::cube E_beta_sq,//J by K by pL1
                                          List anc1, // length pL1; each element is a vector.
                                          arma::vec leaf_desc
){
  // follow tree1:
  int n_desc= leaf_desc.size();
  for (int j=0;j<n_desc;j++){
    int v = leaf_desc(j)-1;
    E_beta.slice(v) *= 0.0;
    E_beta_sq.slice(v) *= 0.0;
    arma::vec curr_anc = anc1[v];
    int n_anc=curr_anc.size();
    for (int u=0;u<n_anc;u++){
      int uu =  curr_anc(u)-1;
      E_beta.slice(v)    += prob1(uu)*mu_gamma.slice(uu);
      E_beta_sq.slice(v) += prob1(uu)*(sigma_gamma.slice(uu)+(1.0-prob1(uu))*pow(mu_gamma.slice(uu),2.0)); //not yet.
    }
    E_beta_sq.slice(v)   += pow(E_beta.slice(v),2.0);
  }
  return List::create(Named("E_beta")=E_beta,
                      Named("E_beta_sq")=E_beta_sq
  );
}

//' Calculate variational moments during the updates
//'
//' Get all moments that need updating when iterating over internal and leaf nodes;
//' for both trees. When updating a node u, all the moments
//' of the descedant nodes will be changed. A recalculation of the moments are necessary
//' when moving on to another node.
//'
//' @param prob2
//' `prob2` is for `s_cu` - a matrix `pL1` by `p2`; in R, a list of pL1 length - each element being of length `p2`.
//' @param mu_alpha variational Gaussian mean vectors (for \code{s_cu=1} component) -
//' this is a pL1 by K-1 by p2 array; in R, we used a list of p2 matrices (each of dimension pL1 by K-1)
//' @param sigma_alpha variational Gaussian variances (for \code{s_cu=1} component)
//' - this is an array of dimension (pL1, K-1, p2); in R, we used a list of p2 matrices,
//' each of dimension pL1 by K-1.
//' @param E_eta,E_eta_sq the intermediate moments to be updated.
//' @param anc2 `anc2` is a list of pL1 vectors, each vector has the node ids of the ancestors in tree2;
//' @param leaf_desc leaf descendants under consideration
//'
//' @return a List
//'
//' \describe{
//'   return List::create(Named("E_eta")=E_eta,
//'    Named("E_eta_sq")=E_eta_sq);
//'}
//'
//' @useDynLib doubletree
//' @importFrom Rcpp sourceCpp
//' @export
// [[Rcpp::export]]
List get_moments_cpp_eco_alpha_doubletree(
    arma::mat prob2,// pL1 by p2
    arma::cube mu_alpha, // pL1 by K-1 by p2
    arma::cube sigma_alpha,//  pL1 by K-1 by p2
    arma::cube E_eta,// pL1,K-1,pL2
    arma::cube E_eta_sq,// pL1,K-1,pL2
    List anc2, // length pL2; each element is a vector.
    arma::vec  leaf_desc
){
  // follow tree2:
  int n_desc= leaf_desc.size();
  for (int j=0;j<n_desc;j++){
    int v = leaf_desc(j)-1;
    E_eta.slice(v) *= 0.0;
    E_eta_sq.slice(v) *= 0.0;
    arma::vec curr_anc = anc2[v];
    int n_anc=curr_anc.size();
    for (int u=0;u<n_anc;u++){
      int uu =  curr_anc(u)-1;
      E_eta.slice(v)        += mtv(mu_alpha.slice(uu),prob2.col(uu));
      E_eta_sq.slice(v)     += mtv(sigma_alpha.slice(uu)+
        mtv(pow(mu_alpha.slice(uu),2.0),(1.0-prob2.col(uu))),
        prob2.col(uu));
    }
    E_eta_sq.slice(v)        += pow(E_eta.slice(v),2.0);
  }
  // return results:
  return List::create(Named("E_eta")=E_eta,
                      Named("E_eta_sq")=E_eta_sq
  );
}
//' Calculate the F term; an array
//'
//' NB: does this work for the case where missing tree1 label can be in multiple domains,
//' and missing/non-missing can happen at the same time for any domain?
//'
//' This function updates the N by pL1 matrix \code{emat} in the package
//' @param psi,g_psi,phi,g_phi local variational parameters
//' @param X transformed data: `2Y-1`; with potential missing responses
//' @param ind_obs_i List of length N; each element is a vector of integers; variable lengths.
//' @param rmat matrix (N by K); row sums are one - the class probabilities for each observation.
//' @param E_beta,E_beta_sq,E_eta,E_eta_sq moment updates produced by \code{\link{get_moments_cpp_doubletree}}
//' @param v1_lookup_NA_replaced pL1+1 for `NA` in tree1 leaf
//' @param v2_lookup `v2_lookup` is a vector of length equal to the total number of rows in \code{X};
//' each element is an integer, indicating which leaf does the observation belong to in tree2
//'
//' @return  n by pL1 by K
//'
//' @useDynLib doubletree
//' @importFrom Rcpp sourceCpp
//' @export
// [[Rcpp::export]]
arma::cube F_doubletree(arma::cube psi, // pL1 by J by K
                        arma::cube g_psi,
                        arma::cube phi, // pL1 by pL2 by K-1
                        arma::cube g_phi,
                        arma::mat X,// N by J; with missing data; <------ do we need to add information about missing index?
                        List ind_obs_i,// length N list; each element is a vector of integers; variable lengths.
                        arma::mat rmat,//N by K,
                        arma::cube E_beta,//J,K,pL1
                        arma::cube E_beta_sq,//J,K,pL1
                        arma::cube E_eta, //pL1,K-1,pL2
                        arma::cube E_eta_sq,//pL1,K-1,pL2
                        arma::vec v1_lookup_NA_replaced,
                        arma::vec v2_lookup){// N by 1; known!
  int n = X.n_rows, J = X.n_cols, K = psi.n_slices, pL1 = psi.n_rows;
  arma::cube res(n,pL1,K);res.zeros();
  for (int i=0;i<n;i++){
    arma::vec curr_ind_obs_i = ind_obs_i[i]; // get the indices of items for which subject i provided a response.
    int J_i = curr_ind_obs_i.size(); // number of observed responses for subject i.
    int v2  = v2_lookup(i)-1;// get leaf id in tree2.
    int v1_NA = v1_lookup_NA_replaced(i)-1;
    if (v1_NA<pL1){//if not missing.
      int v1 = v1_NA; // for the one corresponding to the observed leaf node in tree1.
      for (int k=0;k<K;k++){// for every class k.
        for (int j_obs=0;j_obs<J_i;j_obs++){
          int j = curr_ind_obs_i(j_obs)-1;
          res(i,v1,k) += logexpit_cpp(psi(v1,j,k))+(1.0*X(i,j)*E_beta(j,k,v1)-psi(v1,j,k))*0.5-g_psi(v1,j,k)*(E_beta_sq(j,k,v1)-pow(psi(v1,j,k),2.0));
        }
        if (k<1){//first segment
          res(i,v1,k)   += logexpit_cpp(phi(v1,v2,k))+(E_eta(v1,k,v2)-phi(v1,v2,k))*0.5-g_phi(v1,v2,k)*(E_eta_sq(v1,k,v2)-pow(phi(v1,v2,k),2.0));
        } else if (k< K-1){ // not the first, not the last.
          for (int m=0;m<k-1;m++){
            res(i,v1,k)  += logexpit_cpp(phi(v1,v2,m))+(-E_eta(v1,m,v2)-phi(v1,v2,m))*0.5-g_phi(v1,v2,m)*(E_eta_sq(v1,m,v2)-pow(phi(v1,v2,m),2.0));
          }
          res(i,v1,k)  += logexpit_cpp(phi(v1,v2,k))+(E_eta(v1,k,v2)-phi(v1,v2,k))*0.5-g_phi(v1,v2,k)*(E_eta_sq(v1,k,v2)-pow(phi(v1,v2,k),2.0));
        } else{// k==K-1; the last segment:
          for (int s=0;s<K-1;s++){
            res(i,v1,k)  += logexpit_cpp(phi(v1,v2,s))+(-E_eta(v1,s,v2)-phi(v1,v2,s))*0.5-g_phi(v1,v2,s)*(E_eta_sq(v1,s,v2)-pow(phi(v1,v2,s),2.0));
          }
        }
      }
    }else {// if missing.
      for (int v1=0;v1<pL1;v1++){// for every leaf node in tree1.
        for (int k=0;k<K;k++){// for every class k.
          for (int j_obs=0;j_obs<J_i;j_obs++){
            int j = curr_ind_obs_i(j_obs)-1;
            res(i,v1,k) += logexpit_cpp(psi(v1,j,k))+(1.0*X(i,j)*E_beta(j,k,v1)-psi(v1,j,k))*0.5-g_psi(v1,j,k)*(E_beta_sq(j,k,v1)-pow(psi(v1,j,k),2.0));
          }
          if (k<1){//first segment
            res(i,v1,k)   += logexpit_cpp(phi(v1,v2,k))+(E_eta(v1,k,v2)-phi(v1,v2,k))*0.5-g_phi(v1,v2,k)*(E_eta_sq(v1,k,v2)-pow(phi(v1,v2,k),2.0));
          } else if (k< K-1){ // not the first, not the last.
            for (int m=0;m<k-1;m++){
              res(i,v1,k)  += logexpit_cpp(phi(v1,v2,m))+(-E_eta(v1,m,v2)-phi(v1,v2,m))*0.5-g_phi(v1,v2,m)*(E_eta_sq(v1,m,v2)-pow(phi(v1,v2,m),2.0));
            }
            res(i,v1,k)  += logexpit_cpp(phi(v1,v2,k))+(E_eta(v1,k,v2)-phi(v1,v2,k))*0.5-g_phi(v1,v2,k)*(E_eta_sq(v1,k,v2)-pow(phi(v1,v2,k),2.0));
          } else{// k==K-1; the last segment:
            for (int s=0;s<K-1;s++){
              res(i,v1,k)  += logexpit_cpp(phi(v1,v2,s))+(-E_eta(v1,s,v2)-phi(v1,v2,s))*0.5-g_phi(v1,v2,s)*(E_eta_sq(v1,s,v2)-pow(phi(v1,v2,s),2.0));
            }
          }
        }
      }
    }
  }
  // return results:
  return res;
}

//' Update the variational probabilities of each observation in one of `pL1` leaves in tree1
//'
//' NB: does this work for the case where missing tree1 label can be in multiple domains,
//' and missing/non-missing can happen at the same time for any domain?
//'
//' This function updates the N by pL1 matrix \code{emat} in the package
//' @param curr_F F from [F_doubletree()]
//' @param rmat matrix (N by K); row sums are one - the class probabilities for each observation.
//' @param digamma_emat matrix(pL1 by pL2)
//' @param v2_lookup `v2_lookup` is a vector of length equal to the total number of rows in \code{X};
//' each element is an integer, indicating which leaf does the observation belong to in tree2
//'
//' @return  N by pL1 variational multinomial probabilities for
//' the posterior distribution of N people's membership in the pL1 leaves in tree1;
//' row sums are 1s. `F` is to prepare for `rmat` updates. For subjects with
//' observed leaf label in tree1, the corresponding rows of `res` should be further modified.
//'
//' @useDynLib doubletree
//' @importFrom Rcpp sourceCpp
//' @export
// [[Rcpp::export]]
arma::mat update_emat_with_F_doubletree(arma::cube curr_F, // N by pL1 by K
                                        arma::mat rmat,//N by K
                                        arma::mat digamma_emat,//pL1 by pL2
                                        arma::vec v2_lookup// N by 1; known!
){
  int n = curr_F.n_rows, K = curr_F.n_slices, pL1 = curr_F.n_cols;
  // X = X+0.0;
  arma::mat res(n,pL1);res.zeros();
  arma::vec tmp(n);tmp.zeros();
  arma::mat rF_term(n,pL1);rF_term.zeros();
  for (int i=0;i<n;i++){
    int v2 = v2_lookup(i)-1;
    for (int v=0;v<pL1;v++){
      for (int k=0;k<K;k++){
        rF_term(i,v) += rmat(i,k)*curr_F(i,v,k);//sum over leaf nodes in tree1.
      }
      rF_term(i,v) += digamma_emat(v,v2);
    }
    for (int vv=0;vv<pL1;vv++){
      res(i,vv) = exp(rF_term(i,vv)-logsumexp_row(rF_term.row(i)));
    }
  }
  return res;
}


// //' Update the variational probabilities of each observation in one of `pL1` leaves in tree1
// //'
// //' NB: does this work for the case where missing tree1 label can be in multiple domains,
// //' and missing/non-missing can happen at the same time for any domain?
// //'
// //' This function updates the N by pL1 matrix \code{emat} in the package
// //' @param psi,g_psi,phi,g_phi local variational parameters
// //' @param X transformed data: `2Y-1`; with potential missing responses
// //' @param ind_obs_i List of length N; each element is a vector of integers; variable lengths.
// //' @param rmat matrix (N by K); row sums are one - the class probabilities for each observation.
// //' @param digamma_emat matrix(pL1 by pL2)
// //' @param E_beta,E_beta_sq,E_eta,E_eta_sq moment updates produced by \code{\link{get_moments_cpp_doubletree}}
// //' @param v2_lookup `v2_lookup` is a vector of length equal to the total number of rows in \code{X};
// //' each element is an integer, indicating which leaf does the observation belong to in tree2
// //'
// //' @return  N by pL1 variational multinomial probabilities for
// //' the posterior distribution of N people's membership in the pL1 leaves in tree1;
// //' row sums are 1s. `F` is to prepare for `rmat` updates. For subjects with
// //' observed leaf label in tree1, the corresponding rows of `res` should be further modified.
// //'
// //' @useDynLib doubletree
// //' @importFrom Rcpp sourceCpp
// //' @export
// // [[Rcpp::export]]
// List update_emat_doubletree(arma::cube psi, // pL1 by J by K
//                             arma::cube g_psi,
//                             arma::cube phi, // pL1 by pL2 by K-1
//                             arma::cube g_phi,
//                             arma::mat X,// N by J; with missing data; <------ do we need to add information about missing index?
//                             List ind_obs_i,// length N list; each element is a vector of integers; variable lengths.
//                             arma::mat rmat,//N by K,
//                             arma::mat digamma_emat,//pL1 by pL2
//                             arma::cube E_beta,//J,K,pL1
//                             arma::cube E_beta_sq,//J,K,pL1
//                             arma::cube E_eta, //pL1,K-1,pL2
//                             arma::cube E_eta_sq,//pL1,K-1,pL2
//                             arma::vec v2_lookup){// N by 1; known!
//   int n = X.n_rows, J = X.n_cols, K = psi.n_slices, pL1 = digamma_emat.n_rows;
//   int v = 0;
//   arma::cube F_term(n,pL1,K);F_term.zeros();
//   arma::mat res(n,pL1);res.zeros(); // n is total sample size.
//   arma::vec rF_term(K);rF_term.zeros();
//
//   arma::mat tmp(n,pL1);tmp.zeros();
//   for (int i=0;i<n;i++){
//     arma::vec curr_ind_obs_i = ind_obs_i[i]; // get the indices of items for which subject i provided a response.
//     int J_i = curr_ind_obs_i.size(); // number of observed responses for subject i.
//     int v2  = v2_lookup(i)-1;// get leaf id in tree2.
//     for (int v1=0;v1<pL1;v1++){// for every leaf node in tree1.
//       for (int k=0;k<K;k++){// for every class k.
//         rF_term(k)=0;
//         for (int j_obs=0;j_obs<J_i;j_obs++){
//           int j = curr_ind_obs_i(j_obs)-1;
//           F_term(i,v1,k) += -log(1.0+exp(-psi(v1,j,k)))+(1.0*X(i,j)*E_beta(j,k,v1)-psi(v1,j,k))*0.5-g_psi(v1,j,k)*(E_beta_sq(j,k,v1)-pow(psi(v1,j,k),2.0));
//         }
//         if (k<1){//first segment
//           F_term(i,v1,k)   += -log(1.0+exp(-phi(v1,v2,k)))+(E_eta(v1,k,v2)-phi(v1,v2,k))*0.5-g_phi(v1,v2,k)*(E_eta_sq(v1,k,v2)-pow(phi(v1,v2,k),2.0));
//         } else if (k< K-1){ // not the first, not the last.
//           for (int m=0;m<k-1;m++){
//             F_term(i,v1,k)  += -log(1.0+exp(-phi(v1,v2,m)))+(-E_eta(v1,m,v2)-phi(v1,v2,m))*0.5-g_phi(v1,v2,m)*(E_eta_sq(v1,m,v2)-pow(phi(v1,v2,m),2.0));
//           }
//           F_term(i,v1,k)  += -log(1.0+exp(-phi(v1,v2,k)))+(E_eta(v1,k,v2)-phi(v1,v2,k))*0.5-g_phi(v1,v2,k)*(E_eta_sq(v1,k,v2)-pow(phi(v1,v2,k),2.0));
//         } else{// k==K-1; the last segment:
//           for (int s=0;s<K-1;s++){
//             F_term(i,v1,k)  += -log(1.0+exp(-phi(v1,v2,s)))+(-E_eta(v1,s,v2)-phi(v1,v2,s))*0.5-g_phi(v1,v2,s)*(E_eta_sq(v1,s,v2)-pow(phi(v1,v2,s),2.0));
//           }
//         }
//         rF_term(k) = rmat(i,k)*F_term(i,v1,k);
//       }
//       tmp(i,v1) = digamma_emat(v1,v2)+sum(rF_term);
//     }
//     for (int v1=0;v1<pL1;v1++){
//       res(i,v1) = exp(tmp(i,v1)-logsumexp_row(tmp.row(i)));//we would only need rows corresponding to
//       //subjects with unknown leaf labels in tree1; all 1s if all observed.
//     }
//   }
//   // return results:
//   return List::create(Named("res")=res,
//                       Named("F")=F_term // definitely contains useless terms...for observed tree1 labels.
//   );
// }


//' Update the variational probabilities of each observation in one of `K` classes
//'
//' This function updates the N by K matrix \code{rmat} in the package
//'
//' @param curr_F from \code{\link{F_doubletree}}
//' @param emat matrix (N by K); row sums are one - the class probabilities for each observation.
//'
//' @return  N by K variational multinomial probabilities for
//' the posterior distribution of N people's membership in the K classes in tree1;
//' row sums are 1s.
//'
//' @useDynLib doubletree
//' @importFrom Rcpp sourceCpp
//' @export
// [[Rcpp::export]]
arma::mat update_rmat_with_F_doubletree(arma::cube curr_F,//n,pL1,K
                                        arma::mat emat//N by pL1; in R, should check subjects not in unknown_ids have all 1 emat values.
){
  int n = curr_F.n_rows, K = curr_F.n_slices, pL1 = curr_F.n_cols;
  // X = X+0.0;
  arma::mat res(n,K);res.zeros();
  arma::vec tmp(n);tmp.zeros();
  arma::mat eF_term(n,K);eF_term.zeros();
  for (int i=0;i<n;i++){
    for (int k=0;k<K;k++){
      for (int v=0;v<pL1;v++){
        eF_term(i,k) += emat(i,v)*curr_F(i,v,k);//sum over leaf nodes in tree1.
      }
    }
    for (int k=0;k<K;k++){
      res(i,k) = exp(eF_term(i,k)-logsumexp_row(eF_term.row(i)));
    }
  }
  // return results:
  return res;
}

//' Update gamma's variational distribution. Update the variational mean and variance for logit of
//' class-specific response probabilities (for the \code{s*_u=1} component in tree1.)
//'
//' @param u node id; internal or leaf node in tree1
//' @param g_psi g of local variational parameters
//' @param tau_1_t_u variational Gaussian variances for gamma
//' @param E_beta,E_zeta_u moment updates produced by \code{\link{get_moments_cpp}};
//' \code{E_zeta_u} is directly calculated: \code{prob1[u]*sigma_gamma[u,,]}
//' @param X_zeropad transformed data: `2Y-1`; contains potential missing data.
//' @param rmat a matrix of variational probabilities of all observations
//' belong to K classes; N by K; each row sums to 1
//' @param emat a matrix of variational probability for all observations
//' belonging to pL1 leaf nodes; N by pL1; each row sums to 1. Importantly,
//' for rows with obsered leaf nodes in tree1, we just have an one-hot represention
//' of that cause.
//' @param h_pau a numeric vector of length p indicating the branch length
//' between a node and its parent
//' @param subject_ids_nonmissing the ids of subjects in the leaf descendants of node u;
//' a list of length J, each is a list of subjects nested under u AND have complete info.
//' @param leaf_desc a vector of leaf descendants nested under node `u`
//'
//' @return  a list
//' \describe{
//'   \item{resA}{actually 1/A in the paper, this is variance}
//'   \item{resB}{}
//'   \item{logresBsq_o_A}{}
//' }
//'
//' @useDynLib doubletree
//' @importFrom Rcpp sourceCpp
//' @export
// [[Rcpp::export]]
List update_gamma_subid_doubletree(int u,
                                   arma::cube g_psi,// pL1 by J by K
                                   double tau_1_t_u,
                                   arma::cube E_beta,// J by K by pL1
                                   arma::mat E_zeta_u,//J by K
                                   arma::mat X_zeropad,// N by J
                                   arma::mat rmat,//N by K
                                   arma::mat emat,//N by pL1
                                   arma::vec h_pau,//p1
                                   arma::vec leaf_desc
){
  int n = X_zeropad.n_rows, J = X_zeropad.n_cols, K = g_psi.n_slices;
  int pL1 = g_psi.n_rows;
  int uu = (int) u-1;
  arma::mat resA(J,K);resA.zeros(); // inv A, or variance
  arma::mat resB(J,K);resB.zeros();
  arma::mat logresBsq_o_A(J,K);logresBsq_o_A.zeros();

  X_zeropad = 1.0*X_zeropad;

  int n_leaf_desc = leaf_desc.size();
  for (int j=0;j<J;j++){
    for (int k=0;k<K;k++){
      resA(j,k) = 1.0/(tau_1_t_u*h_pau(uu));
      for (int v1=0;v1<n_leaf_desc;v1++){
        int curr_v  = leaf_desc(v1)-1;
        resA(j,k)  += 2.0*g_psi(curr_v,j,k)*sum(emat.col(curr_v)%rmat.col(k));
        resB(j,k)  += sum(emat.col(curr_v)%(rmat.col(k)%(X_zeropad.col(j)*0.5-
          X_zeropad.col(j)%X_zeropad.col(j)*2.0*g_psi(curr_v,j,k)*(E_beta(j,k,curr_v)-E_zeta_u(j,k)))));
      }
      logresBsq_o_A(j,k) = 2.0*log(abs(resB(j,k)))-log(resA(j,k));
      resA(j,k) = 1.0/resA(j,k);
    }
  }
  return List::create(Named("resA")=resA,
                      Named("resB")=resB,
                      Named("logresBsq_o_A")=logresBsq_o_A);
}



//get_moments_cpp_eco

//' Update alpha's variational distribution.
//'
//' @param u node id; internal or leaf node in tree1
//' @param v1 leaf node id in tree1.
//' @param g_psi,g_phi g of local variational parameters
//' @param tau_2_t_u variational Gaussian variances for gamma
//' @param E_eta,E_xi_u moment updates produced by \code{\link{get_moments_cpp}};
//' \code{E_xi_u} is directly calculated
//' @param X transformed data: `2Y-1`; contains potential missing data.
//' @param rmat a matrix of variational probabilities of all observations
//' belong to K classes; N by K; each row sums to 1
//' @param emat a matrix of variational probability for all observations
//' belonging to pL1 leaf nodes; N by pL1; each row sums to 1. Importantly,
//' for rows with obsered leaf nodes in tree1, we just have an one-hot represention
//' of that cause.
//' @param h_pau a numeric vector of length p indicating the branch length
//' between a node and its parent; for tree2
//' @param levels a vector of possibly repeating integers from 1 to Fg2
//' @param subject_ids integer ids for subjects nested under node `u`
//' @param v2_lookup a vector of length equal to the total number of rows in X;
//' each element is an integer, indicating which leaf does the observation belong to in tree2.
//'
//' @return  a list
//' \describe{
//'   \item{resC}{actually 1/C in the paper, this is variance}
//'   \item{resD}{}
//'   \item{logresDsq_o_C}{}
//' }
//'
//' @useDynLib doubletree
//' @importFrom Rcpp sourceCpp
//' @export
// [[Rcpp::export]]
List update_alpha_subid_doubletree(
    int u,
    int v1,
    arma::cube g_phi,// pL1 by pL2 by K-1
    double tau_2_t_u,
    arma::cube E_eta,//pL1 by K-1 by pL2
    arma::mat E_xi_u,// pL1 by K-1
    arma::mat X,// N by J
    arma::mat rmat,//N by K
    arma::mat emat,//N by pL1
    arma::vec h_pau,//p2
    arma::vec levels,//p2
    arma::vec subject_ids,
    arma::vec v2_lookup//N by 1
){
  int n = subject_ids.size(), J = X.n_cols, K = rmat.n_cols;
  int p = h_pau.size(), pL1 = E_eta.n_rows;
  int ii = 0;
  int vv = 0;
  int uu = (int) u-1;

  arma::vec resC(K-1);resC.zeros(); // inv C, or variance
  arma::vec resD(K-1);resD.zeros();
  arma::vec pre_resC(2);pre_resC.zeros();
  arma::vec logresDsq_o_C(K-1);logresDsq_o_C.zeros();

  int vv1 = v1-1;
  for (int k=0;k<K-1;k++){
    resC(k) = -1.0*log(tau_2_t_u)-1.0*log(h_pau(uu));
    // resC(k) = 1/(tau_2_t_u*h_pau(uu));
    for (int i=0;i<n;i++){
      ii = (int) subject_ids(i)-1;
      vv = (int) v2_lookup(ii)-1;
      for (int m=k;m<K;m++){
        pre_resC(0) = resC(k);
        pre_resC(1) = log(2.0)+log(rmat(ii,m))+log(g_phi(vv1,vv,k))+log(emat(ii,vv1));
        resC(k) = logsumexp(pre_resC);// need to make sure emat have correct one-hot representations for observed CODs.
        // resC(k) += 2.0*rmat(ii,m)*g_phi(vv1,vv,k)*emat(ii,vv1);// need to make sure emat have correct one-hot representations for observed CODs.
        if (m<k+1){
          resD(k) += emat(ii,vv1)*(rmat(ii,m)*0.5 - 2.0*rmat(ii,m)*g_phi(vv1,vv,k)*(E_eta(vv1,k,vv)-E_xi_u(vv1,k)));
        }else{
          resD(k) += emat(ii,vv1)*(-rmat(ii,m)*0.5 - 2.0*rmat(ii,m)*g_phi(vv1,vv,k)*(E_eta(vv1,k,vv)-E_xi_u(vv1,k)));
        }
      }
    }
    logresDsq_o_C(k) = 2.0*log(abs(resD(k)))-resC(k);
    resC(k) = expm1(-resC(k))+1.0;
  }
  // resC = 1.0/resC;
  return List::create(Named("resC")=resC,// this is actually 1/C in the paper, the variance
                      Named("resD")=resD,
                      Named("logresDsq_o_C")=logresDsq_o_C);
}


// others not organized:
// get_est_cpp: summary
// some functions to calculates the ELBO components.

//' calculate line 1, 2, and 15 of ELBO*
//'
//' Intended to be faster than its counterpart implementation in `R`.
//'
//' This function updates the N by pL1 matrix \code{emat} in the package
//' @param F intermediate values from updates of rmat and emat; dimension: N by pL1 by K;
//' currently for observations with known leaf labels in tree1 AND tree2, there is redundancy
//' because once we know i, we know the leaf 1 label, hence no need for
//' the second dimension
//' @param digamma_emat matrix(pL1 by pL2)
//' @param rmat matrix (N by K); row sums are one - the class probabilities for each observation.
//' @param emat matrix (N by pL1)
//' @param v1_lookup_NA_replaced vector; includes values of pL+1 for unknown leaf label in tree1.
//' @param v2_lookup `v2_lookup` is a vector of length equal to the total number of rows in \code{X};
//' each element is an integer, indicating which leaf does the observation belong to in tree2
//'
//' @return  line 1, 2, 15 of the ELBO* in Appendix
//'
//' @useDynLib doubletree
//' @importFrom Rcpp sourceCpp
//' @export
// [[Rcpp::export]]
List get_line1_2_15_doubletree(arma::cube F,//N by pL1 by K,
                               arma::mat digamma_emat,//pL1 by pL2,
                               arma::mat rmat,
                               arma::mat emat,
                               arma::vec v1_lookup_NA_replaced,// N by 1; includes values of pL+1 for unknown leaf id in tree1.
                               arma::vec v2_lookup){// N by 1; known!
  int n = F.n_rows, K = F.n_slices, pL1 = F.n_cols;
  double res1 = 0.0;
  double res2 = 0.0;
  double res3 = 0.0;

  for (int i=0;i<n;i++){
    int v1 = v1_lookup_NA_replaced(i)-1;//get leaf id in tree1; pL1+1 (pL in C indexing) represents missing leaf 1 label.
    int v2 = v2_lookup(i)-1;// get leaf id in tree2.
    if (v1<pL1){// not missing
      res1 += digamma_emat(v1,v2);
      for (int k=0;k<K;k++){
        res1 += rmat(i,k)*F(i,v1,k);
        res3 += - xlogx(rmat(i,k));
      }
    }else{
      for (int vv=0;vv<pL1;vv++){
        res3 += -xlogx(emat(i,vv));
        res2 += emat(i,vv)*digamma_emat(vv,v2);
        for (int k=0;k<K;k++){
          res2 += emat(i,vv)*rmat(i,k)*F(i,vv,k);
        }
      }
      for (int k=0;k<K;k++){
        res3 += -xlogx(rmat(i,k));
      }
    }
  }
  // return results:
  return List::create(Named("res1")=res1,
                      Named("res2")=res2,
                      Named("res3")=res3);
}













