########################
# helper functions:
#######################

#' multi-expit (numerical overflow dealt with)
#'
#' @param v a vector of length K
#'
#' @return a vector of probabilities of length K; sums to one
#' @examples
#' sigmoid <- function(v) mexpit(c(0,v))[2]
#' mexpit(c(1,0,0))
#' @export
mexpit  <- function(v) exp(v - matrixStats::logSumExp(v))
#inverse of mexpit: alrInv()


#' sigmoid (copied from moretrees:)
#'
#' @param p a probability
#'
#' @export
#' @family utility function
logit   <- function(p) log(p) - log(1 - p)


#' log(1+exp(x)) (copied from moretrees:)
#'
#' @param x a number
#'
#' @export
#' @family utility function
log1p_exp <- function(x) {
  # computes log(1 + exp(x))
  if (x > 20) {
    return(x)
  } else {
    return(log1p(exp(x)))
  }
}

#' log expit vectorized
#'
#' @param x a vector of numbers
#'
#' @export
#' @family utility function
log1p.exp.vec <- function(x){
  # computes log(1 + exp(x)) for x a vector
  y <- x
  which.small <- x <= 20
  y[which.small] <- log1p(exp(x[which.small]))
  return(y)
}


#' log expit
#'
#' @param x a number or a vector of numbers
#'a
#' @export
#' @family utility function
logexpit <- function(x) {
  # computes -log(1+exp(-x)) for x a vector
  -log1p.exp.vec(-x)
}

#' expit
#'
#' @param x a number
#'
#' @export
#' @family utility function
expit <- function(x) {
  # computes 1/(1+exp(-x)) for x a vector
  exp(logexpit(x))
}


#' g
#'
#' @param eta local variational parameter
#'
#' @export
#' @family VI functions
g_fun <- function(eta){

  ## Inputs ##

  # eta = a variational parameter

  ## Outputs ##

  # g(eta), a function of eta

  ## Code ##
  (1/(2*eta))*(1/(1+exp(-eta))-0.5)
}

#' g 0
#'
#' @param eta local variational parameter
#'
#' @export
#' @family VI functions
g_fun0 <- function(eta){

  ## Inputs ##

  # eta = a variational parameter

  ## Outputs ##

  # g(eta), a function of eta

  ## Code ##
  if (eta ==0){return(1/8)}
  (1/(2*eta))*(1/(1+exp(-eta))-0.5)
}



#' g, vectorized
#'
#' @param eta local variational parameter
#'
#' @export
#' @family VI functions
g_fun.vec <- function(eta){

  ## Inputs ##

  # eta = a numeric vector

  ## Outputs ##

  # g(eta), a numeric vector containing the values of the function g evaluated at each element of eta

  ## Code ##

  D <- dim(as.array(eta))
  y <- array(1/8,D)
  which.nonzero <- eta != 0
  y[which.nonzero] <- g_fun(eta[which.nonzero])
  y
}

#' transform a vector of probabilities that sum to one to stick proportions.
#'
#' @param x a vector of probabilities (K); sums to `1`
#' @return a vector K, with last element of 1; the elements are stick lengths in
#' the remaining part
#'
#' @examples
#'
#' prob2stick(c(0.5,0.2,0.3))
#'
#' @export
prob2stick <- function(x){
  res <- x
  res[1] <- x[1]
  for (i in 2:length(x)){
    res[i] <- x[i]/(1-sum(x[1:(i-1)]))
  }
  res
}


#' create a list with members being the matrices along a specified dimension of an array
#'
#' This is from a StackOverflow answer
#'
#' @param a an array
#' @param n along which dimension to create a list
#'
#' @return a list
#'
#' @examples
#'
#' myarray <- array(c(1,2,3,4,5,6,7,8),c(2,2,2))
#'
#' split_along_dim(myarray,1)
#' split_along_dim(myarray,2)
#' split_along_dim(myarray,3)
#' @importFrom stats setNames
#' @export
#' @family utility functions
split_along_dim <- function(a, n){
  setNames(lapply(split(a, arrayInd(seq_along(a), dim(a))[, n]),
                  array, dim = dim(a)[-n], dimnames(a)[-n]),
           dimnames(a)[[n]])
}


#' lower bound of expit(x)
#'
#' lower bound is quadratic in the exponent
#'
#' @param xi local variational parameter, positive value, this is where
#' the `expit(xi) = lower_bd(xi)`; this is where the quadratic
#' and expit curve contact
#'
#' @return a positive value between `0` and `1`
#'
#' @examples
#'
#' par(mfrow=c(2,3))
#' xi_seq <- c(0.1,0.5,1,2,4,6)
#'   x <- seq(-5,5,by=0.1)
#' for (i in 1:length(xi_seq)){
#'   xi <- xi_seq[i]
#'   y1 <- lower_bd(xi)(x)
#'   y2 <- lower_bd(-xi)(x)
#'
#'   plot(x,y1,type="l",col="red",ylim=c(0,1),main=paste0("xi= ",xi))
#'   points(x,y2,type="l",col="blue")
#'   points(x,1/(1+exp(-x)),type="l",col="orange")
#'   abline(v=c(-xi,xi),col="gray")
#'   legend("topleft",c("expit","bound"),col=c("orange","blue"),lty=c(1,1))
#' }
#'
#' @references
#' \itemize{
#' \item Jaakkola, Tommi S., and Michael I. Jordan. "Bayesian parameter estimation via variational methods." Statistics and Computing 10.1 (2000): 25-37.
#' <http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.399.9368&rep=rep1&type=pdf>
#'}
#' @export
#' @family VI functions
lower_bd <- function(xi){
  function(z){
    expit(xi)*exp((z-xi)/2-g_fun0(xi)*(z^2-xi^2))
  }
}



#' lower bound of a vector of probabilities that sum to one
#'
#' the lower bound is based on the best set of local variational parameters
#' which comprise of logit of the stick-breaking form of the supplied vector
#'
#' @param x a vector of probabilities that sum to one
#'
#' @return approximation to a vector of probabilities
#'
#' @examples
#'
#' # based on Tsiatis 2016 NeuroIPS
#' approx <- function(x){
#'   res = rep(NA,length=length(x))
#'   for (i in 1:length(res)){
#'     curr_v <- x[i] - x[-i]
#'     res[i]  = prod(expit(curr_v))
#'   }
#'   res
#' }
#'
#' tau = rep(0.25,4)
#' barplot(rbind(tau,approx_sb(tau),approx(tau)),beside=TRUE,
#'         legend.text=c("truth","sb + quad (lotR)","1 vs other + quad"),
#'         main="truth 1")
#'
#' tau = c(0.5,0.3,0.15,0.05)
#' barplot(rbind(tau,approx_sb(tau),approx(tau)),beside=TRUE,
#'         legend.text=c("truth","sb + quad (lotR)","1 vs other + quad"),
#'         main="truth 2")
#'
#' @references
#' \itemize{
#' \item Jaakkola, Tommi S., and Michael I. Jordan. "Bayesian parameter estimation via variational methods." Statistics and Computing 10.1 (2000): 25-37.
#' <http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.399.9368&rep=rep1&type=pdf>
#' \item Titsias M(2016). One-vs-each approximation to softmax for scalable estimation of probabilities. Advances in Neural Information Processing Systems.
#' <https://papers.nips.cc/paper/6468-one-vs-each-approximation-to-softmax-for-scalable-estimation-of-probabilities.pdf>
#'}
#' @export
#' @family VI functions
approx_sb <- function(x){
  logit_stick_lengths <-  logit(prob2stick(x))[-length(x)] # best local variational parameters, e.g., psi, phi in lotR.
  res = rep(1,length=length(x))
  hm = rep(1,length=length(x))
  hp = rep(1,length=length(x))
  hp[1] = lower_bd(logit_stick_lengths[1])(logit_stick_lengths[1])
  res[1] = hp[1]
  for (i in 2:(length(res)-1)){
    hm[i] = lower_bd(logit_stick_lengths[i-1])(-logit_stick_lengths[i-1])
    hp[i] = lower_bd(logit_stick_lengths[i])(logit_stick_lengths[i])
    res[i] = prod(hm[1:(i)])*hp[i]
  }
  res[length(x)] = prod(hm[1:(length(x)-1)])*lower_bd(logit_stick_lengths[length(x)-1])(-logit_stick_lengths[length(x)-1])
  res
}


#' Takes any number of R objects as arguments and returns a list whose names are
#' derived from the names of the R objects.
#'
#' Roger Peng's list labeling challenge from
#' <http://simplystatistics.tumblr.com/post/11988685443/computing-on-the-language>.
#' Code copied from <https://gist.github.com/ajdamico/1329117/0134148987859856fcecbe4446cfd37e500e4272>
#'
#' @param ... any R objects
#'
#' @return a list as described above
#'
#' @examples
#' #create three example variables for a list
#' x <- 1
#' y <- 2
#' z <- "hello"
#' #display the results
#' make_list( x , y , z )
#' @export
#' @family utility functions
#create the function
make_list <- function(...) {
  #put all values into a list
  argument_values <- list(...)

  #save all argument names into another list
  argument_names <- as.list(sys.call())

  #cycle through the first list and label with the second, ignoring the function itself
  for (i in 2:length(argument_names)) {
    names(argument_values)[i - 1] <- argument_names[i]
  }

  #return the newly-labeled function
  argument_values
}


# n = 10000
# J = 3
# se = sqrt(9/4)
# mu = 0
# xmat <- matrix(rnorm(n*J,mu,se),nrow=n,ncol=J)
# y <- apply(cbind(xmat,0),1,mexpit)
#
# par(mfrow=c(1,J+1))
# for (s in 1:(J+1)){
#   hist(y,breaks="Scott")
# }
#
# rowMeans(y)



#' generate stick-breaking prior (truncated) from a vector of random probabilities
#'
#' @param u a vector of probabilities, with the last element 1.
#' Each element means the fraction of what is left to be broken. The last
#' element is 1 because we truncate the length of the stick to be `length(u)`
#'
#' @return a vector of the same length as u; sum to 1.
#'
#' @examples
#'
#' graphics::par(mfrow=c(3,3),oma=c(0,1,5,0),
#'    mar=c(1,2,1,1))
#' for (iter in 1:9){
#'  u   <- c(rbeta(9,1,0.8),1)
#'  res <- tsb(u)
#'  barplot(res,ylim=c(0,1),main=paste0("Random Sample #", iter),ylab="Probability")
#' }
#' graphics::mtext("Truncated Stick-Breaking Dist. (10 segments)",3,
#'      outer=TRUE,cex=1.5,line=1.5)
#' @export
#'
tsb <- function(u){
  K <- length(u)
  small_tol <- 1e-200
  if (abs(u[K]-1)>1e-6) {stop("==The last element of u must be 1 for truncated stick-breaking!==\n")}
  w <- rep(NA,K)
  w[1] <- u[1]
  for (k in 2:(K)){
    w[k] <- exp(log(w[k-1]+small_tol)-log(u[k-1]+small_tol)+log(1-(u[k-1]+small_tol))+log(u[k]+small_tol))
  }
  w
}

#' generate random vectors from Dirichlet distribution
#'
#' @param n sample size
#' @param v a vector of positive values
#'
#' @return matrix of n rows and `length(v)` columns; each row sums to 1
#'
#' @examples
#'
#' myrdirich(10,c(1,1,1))
#'
#' @importFrom stats rgamma
#' @export
myrdirich <- function(n,v){
  res <- matrix(NA,nrow=n,ncol=length(v))
  for (j in 1:length(v)){
    res[,j] <- rgamma(n,v[j])
  }
  sweep(res,MARGIN=1,rowSums(res),"/")
}





## 1. top k cause classification accuracy


#' get top k ids for each row
#'
#' @param probs a matrix where columns are causes; the entries are probabilities of causes
#' @param s top s
#'
#' @return data frame; for each row in probs, indicating the top s causes
#'
#' @examples
#'
#' xx <- matrix(c(0.3,0.1,0.2,0.4,
#'                0.4,0.1,0.3,0.2),nrow=2,ncol=4,byrow=TRUE)
#' get_topk_COD(xx,2)
#'
#' @export
#'
get_topk_COD <- function(probs,s){
  # which.k <- function(x, k){which(x == sort(x, decreasing = TRUE)[k])}
  # v <- ncol(probs)
  # res <- data.frame(cause1=apply(probs, 1, which.k, 1))
  # if (s>1){
  #   for (s2 in 2:s){
  #     res <- cbind(res,apply(probs, 1, which.k, s2))
  #   }
  # }
  # names(res) <- paste("cause",1:s,sep="")
  # res
  res <- data.frame(matrix(t(apply(probs,1,function(v) order(v,decreasing = TRUE)[1:s])),ncol=s))
  names(res) <- paste("cause",1:s,sep="")
  res
}

#' get top k accuracy
#'
#' @param pred_topk a matrix of integers; each row indicates top 1, 2, ..., k
#' causes
#' @param truth a vector of integers; its length equals the number of rows in
#' `pred_topk`.
#'
#' @return numeric; between 0 and 1; the higher the more accurate.
#'
#' @examples
#'
#' xx <- matrix(c(0.3,0.1,0.2,0.4,
#'                0.4,0.1,0.3,0.2),nrow=2,ncol=4,byrow=TRUE)
#' pred_top2 <- get_topk_COD(xx,2)
#' acc_topk(pred_top2,c(4,2))
#'
#' @export
acc_topk <- function(pred_topk,truth){
  # res = 0
  # for (k in 1:ncol(pred_topk)){
  #   res <- res + sum(pred_topk[,k] == truth,na.rm=TRUE) / length(truth)
  # }
  # res
  sum(sapply(1:length(truth),function(i) 0+truth[i]%in%pred_topk[i,]))/length(truth)
}

#' heatmap
#'
#' @param mat a matrix
#'
#' @return image
#'
#' @importFrom stats heatmap
#' @export
myhm <- function(mat){
  heatmap(mat,Rowv=NA,Colv=NA)
}

#' sum over signed numbers, with some logsumexp applied to positive and negative
#' values separately
#'
#' @param logabsv a vector; positive or negative; log of abs(v)
#' @param sign sign of v, a vector; +1 or -1.
#' @return image
#'
#' @examples
#' set.seed(123)
#' x <- rnorm(n = 100, mean = 1000, sd = 10)
#' comp <- c(log(sum(exp(x))),signlogsumexp(x,rep(1,length(x)))$res,logsumexp(x))
#' print(comp)
#'
#' sign_indicators <- c(rep(1,length(x)-50),rep(-1,length(x)-50))
#' comp2 <- c(log(sum(sign_indicators*exp(x))),signlogsumexp(x,sign_indicators)$res)
#' print(comp2)
#'
#' x = c(.5,-.4,-.1)
#' sum(x) # not zero
#'
#' signlogsumexp(c(log(abs(x))),sign(x)) # not zero, but much closer to zero.
#'
#' @importFrom matrixStats logSumExp
#'
#' @export
signlogsumexp <- function(logabsv,sign){
  res1 <- NULL
  res2 <- NULL
  if (sum(sign>0)>0){
    logabsv_p <- logabsv[sign>0]
    res1 <- logSumExp(logabsv_p)
  }
  if (sum(sign<0)>0){
    logabsv_m <- logabsv[sign<0]
    res2 <- logSumExp(logabsv_m)
  }

  if (is.null(res1)){
    res      <- logSumExp(logabsv)
    res_sign <- -1
    return(make_list(res,res_sign))
  }
  if (is.null(res2)){
    res      <- logSumExp(logabsv)
    res_sign <- 1
    return(make_list(res,res_sign))
  }
  res_sign <- (res1<res2)*(-2)+1
  a <- max(res1,res2)
  b <- min(res1,res2)
  res <- b+log(expm1(a-b))
  return(make_list(res,res_sign))
}


#' get CSMF accuracy
#'
#' @param csmf a vector that sums to one; positive values
#' @param truth a vector that sums to one; positive values
#' @param undet `NULL` default.
#'
#' @return numeric; between 0 and 1; the higher the more accurate.
#'
#' @export
getCSMFacc <- function (csmf, truth, undet = NULL)
{
  if (methods::is(csmf, "insilico")) {
    if (!is.null(names(truth))) {
      order <- match(colnames(csmf$csmf), names(truth))
      if (is.na(sum(order))) {
        stop("Names not matching")
      }
      truth <- truth[order]
    }
    acc <- 1 - apply(abs(truth - t(csmf$csmf)), 2, sum)/2/(1 -
                                                             min(truth))
  }
  else {
    if (!is.null(undet)) {
      if (undet %in% names(csmf)) {
        csmf <- csmf[-which(names(csmf) == undet)]
      }
      else {
        warning("The undetermined category does not exist in input CSMF.")
      }
    }
    if (!is.null(names(csmf)) & !is.null(names(truth))) {
      order <- match(names(csmf), names(truth))
      if (is.na(sum(order))) {
        stop("Names not matching")
      }
      truth <- truth[order]
    }
    acc <- 1 - sum(abs(truth - csmf))/2/(1 - min(truth))
  }
  return(acc)
}

## 3. RMSE
## 4. aRI for assessing cluster estimation accuracy

