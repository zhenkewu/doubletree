#' update variational parameters
#'
#' Used by [fit_nlcm_doubletree()]; followed by [update_hyperparams_doubletree()]
#'
#' @param leaf_ids_units,leaf_ids_nodes,ancestors,h_pau,levels,v_units,subject_id_list,ind_obs_i,ind_obs_j,scenario
#' outputs from [design_doubletree()]. Among them,
#' `subject_id_list` is currently not used in the initialization
#' function [initialize_nlcm_doubletree()], but are useful for updating functions here.
#' @param X,n,J,p1,p2,pL1,pL2,Fg1,Fg2,cardanc1,cardanc2 data and design information computed from the outputs
#' of [initialize_nlcm_doubletree()]; we use `R.utils::do.Call` to get their values from `dsgn`.
#' @param  prob1,prob2,mu_gamma,mu_alpha,rmat,emat,dirich_mat,sigma_gamma,sigma_alpha,tau_1_t,tau_2_t,a1_t,b1_t,a2_t,b2_t
#' variational parameters
#' @param psi,g_psi,phi,g_phi,tau_1,tau_2 parameters updated by `update_hyperparams_doubletree()`
#' @param a1,b1,a2,b2,dmat,K,LD,s1_u_zeroset,s1_u_oneset,s2_cu_zeroset,s2_cu_oneset
#' @param do_tree1_update `TRUE` to update `mu_gamma`,`sigma_gamma`, `prob1`,`tau_1`
#' fixed hyperparameters
#' @importFrom matrixStats logSumExp
#' @family Internal VI functions
update_vi_params_doubletree <- function(
  leaf_ids_units,
  leaf_ids_nodes,
  ancestors,
  h_pau,
  levels,
  v_units,subject_id_list,ind_obs_i,ind_obs_j,scenario,
  X,n,J,p1,p2,pL1,pL2,Fg1,Fg2,cardanc1,cardanc2,# data and design - outputs or derivatives from design outputs.
  prob1,prob2,mu_gamma,mu_alpha,rmat,emat,dirich_mat,
  sigma_gamma,sigma_alpha,tau_1_t,tau_2_t,a1_t,b1_t,a2_t,b2_t, # vi parameters.
  psi,g_psi,phi,g_phi,tau_1,tau_2, # local variational parameters and prior variances for the slab components.
  a1,b1,a2,b2,dmat,K,LD,s1_u_zeroset,s1_u_oneset,s2_cu_zeroset,s2_cu_oneset, #fixed hyper-parameters; updated.
  do_tree1_update
){
  if (is.null(LD) & K>1){LD <- TRUE}
  X <- as.matrix(X)
  X_zeropad <- X
  X_zeropad[is.na(X)] <- 0
  # calculate some useful quantities:
  v1_units_NA_replaced <- v_units[[1]]
  v1_units_NA_replaced[is.na(v_units[[1]])] <- pL1+1 # this is replacing NA; convenient in C functions.

  if (do_tree1_update){
    ## tree1 restrictions:
    seq_update1 <- 1:p1
    if (!is.null(s1_u_zeroset)){ # not updating nodes that are set to zeros.
      if (length(setdiff(s1_u_zeroset,1:p1))!=0){stop("[doubletree] 's1_u_zeroset has elements not between 1 and p1.'")}
      seq_update1 <- (1:p1)[-s1_u_zeroset]
      prob1[s1_u_zeroset] <- 0
      if (do_tree1_update){
        for (v in s1_u_zeroset){ # this may not matter; but just to be logically clear.
          mu_gamma[[v]] <- mu_gamma[[v]]*0.0               # no diffusion
          sigma_gamma[[v]] <- matrix(tau_1[levels[[1]][v]]*h_pau[[1]][v],nrow=J,ncol=K) # just prior variance.
        }
      }
    }
    if (!is.null(s1_u_oneset)){
      prob1[s1_u_oneset] <- 1
      if (length(setdiff(s1_u_oneset,1:p1))!=0){stop("[doubletree] 's1_u_oneset has elements not between 1 and p1.'")}
    }
  }

  ## tree2 restrictions:
  seq_update2 <- vector(mode = "list", length = pL1)
  for (l in 1:pL1){seq_update2[[l]] <- 1:p2}
  if (!is.null(s2_cu_zeroset)){ # not updating nodes that are set to zeros.
    for (v in 1:pL1){
      if (!is.null(s2_cu_zeroset[[v]])){
        if (length(setdiff(s2_cu_zeroset[[v]],1:p2))!=0){
          stop(paste0("[doubletree] 's2_cu_zeroset has elements not between 1 and p2; check c= '",v,"."))}
        seq_update2[[v]] <- (1:p2)[-s2_cu_zeroset[[v]]]
        prob2[[v]][s2_cu_zeroset[[v]]] <- 0
        for (u in s2_cu_zeroset[[v]]){ # this may not matter; but just to be logically clear.
          mu_alpha[[u]] <- mu_alpha[[u]]*0.0                                     # no diffusion.
          sigma_alpha[[u]] <- matrix(tau_2[levels[[2]][u]]*h_pau[[2]][u],nrow=pL1,ncol=K-1) # just prior variance.
        }
      }
    }
  }
  if (!is.null(s2_cu_oneset)){
    for (v in 1:pL1){
      if (!is.null(s2_cu_zeroset[[v]])){
        prob2[[v]][s2_cu_oneset[[v]]] <- 1
        if (length(setdiff(s2_cu_oneset[[v]],1:p2))!=0){
          stop(paste0("[doubletree] 's2_cu_oneset has elements not between 1 and p2; check c= '",v,"."))}
      }
    }
  }

  if (!exists("E_beta_sq") || !exists("E_eta_sq") || !exists("E_beta") || !exists("E_eta")){
    # calculate initial moments that are required in the VI updates (do so only when not available):
    moments_cpp <- get_moments_cpp_doubletree(
      prob1,as.matrix(do.call("rbind",prob2)),
      array(unlist(mu_gamma),c(J,K,p1)),array(unlist(sigma_gamma),c(J,K,p1)),
      array(unlist(mu_alpha),c(pL1,K-1,p2)),array(unlist(sigma_alpha),c(pL1,K-1,p2)),
      ancestors[[1]],ancestors[[2]],
      cardanc1,cardanc2)
    # needed in updating rmat, and for update_hyperparams:
    E_beta_sq      <- moments_cpp$E_beta_sq
    E_eta_sq       <- moments_cpp$E_eta_sq
    E_beta         <- moments_cpp$E_beta
    E_eta          <- moments_cpp$E_eta
  }

  # update mu_gamma,sigma_gamma,prob1:----------------------------------------
  if (do_tree1_update){
    tau_1_t <- tau_1[levels[[1]]]
    print("tree1_updated.")
    for (u1 in seq_update1){
      if (!is.null(s1_u_oneset) && u1%in%s1_u_oneset){
        prob1[u1] <- 1
      }
      gamma_update <- update_gamma_subid_doubletree(
        u1,g_psi,
        tau_1_t[u1],
        E_beta,as.matrix(prob1[u1]*mu_gamma[[u1]],nrow=J,ncol=K),
        as.matrix(X_zeropad),
        rmat,emat,
        h_pau[[1]],
        leaf_ids_nodes[[1]][[u1]])

      mu_gamma[[u1]]    <-   gamma_update$resB*gamma_update$resA #  J by K
      sigma_gamma[[u1]] <-   gamma_update$resA

      # update probability:---------------------------------------------------
      w1_u <- digamma(a1_t[levels[[1]][u1]])-digamma(b1_t[levels[[1]][u1]])+
        0.5*sum(gamma_update$resBsq_o_A)-
        0.5*J*K*log(tau_1_t[u1]*h_pau[[1]][u1])-0.5*sum(-log(sigma_gamma[[u1]]))

      prob1[u1] <- expit(w1_u)
      if (!is.null(s1_u_oneset) && u1%in%s1_u_oneset){
        prob1[u1] <- 1
      }
      # recalculate moments; the ones that matters are the descedant nodes of u1:
      moments_cpp <- get_moments_cpp_eco_gamma_doubletree(
        prob1,array(unlist(mu_gamma),c(J,K,p1)),array(unlist(sigma_gamma),c(J,K,p1)),
        E_beta,E_beta_sq,ancestors[[1]],leaf_ids_nodes[[1]][[u1]])
      E_beta <- moments_cpp$E_beta
      E_beta_sq <- moments_cpp$E_beta_sq
    }
  }

  # update for the mu_alpha, sigma_alpha, prob2--------------------------------
  tau_2_t <- tau_2[levels[[2]]]
  for (v1 in 1:pL1){
    if (!LD){prob2[[v1]]<- c(1,rep(0,p2-1))} # no updates for prob2 if conditional independence is assumed.
    for (u2 in seq_update2[[v1]]){
      if (LD){ # only update alpha related VI parameters when K=2, and LD being TRUE.
        # update mu_alpha:---------------------------------------------------
        if (!is.null(s2_cu_oneset[[v1]]) && u2%in%s2_cu_oneset[[v1]]){
          prob2[[v1]][u2] <- 1
        }
        # alpha_update <- update_alpha_subid_doubletree(
        #   u2,v1,g_phi,tau_2_t[u2],E_eta,
        #   as.matrix(sweep(mu_alpha[[u2]],MARGIN=1,do.call("rbind",prob2)[,u2],"*")),
        #   as.matrix(X),rmat,emat,h_pau[[2]],levels[[2]],
        #   subject_id_list[[2]][[u2]],v_units[[2]])

        alpha_update <- update_alpha_subid_doubletree0(
          u2,v1,g_phi,tau_2_t[u2],E_eta,
          as.matrix(sweep(mu_alpha[[u2]],MARGIN=1,do.call("rbind",prob2)[,u2],"*")),
          as.matrix(X),rmat,emat,h_pau[[2]],levels[[2]],
          subject_id_list[[2]][[u2]],v_units[[2]])

        # tmp_resD <- alpha_update$resD_mat
        # absv <- abs(tmp_resD)
        # signv <- sign(tmp_resD)

        # curr_resD <- rep(NA,K-1)
        # for (k in 1:(K-1)){
        #   tmpres <- signlogsumexp(log(absv[signv[,k]!=0,k]),
        #                           signv[signv[,k]!=0,k])
        #   curr_resD[k] <- tmpres$res_sign*exp(tmpres$res)
        # }

        #mu_alpha[[u2]][v1,] <- curr_resD*alpha_update$resC
        mu_alpha[[u2]][v1,]    <-  alpha_update$resD*alpha_update$resC
        sigma_alpha[[u2]][v1,] <-  alpha_update$resC

        # print(sum(curr_resD-alpha_update0$resD))
        # print(sum(alpha_update$resC-alpha_update0$resC))


        # update prob2 ---------------------------------------------------
        w2_cu <- digamma(a2_t[v1,levels[[2]][u2]])-digamma(b2_t[v1,levels[[2]][u2]])+
          0.5*sum(alpha_update$resDsq_o_C)-
          # 0.5*sum(curr_resD^2*alpha_update$resC)-
          0.5*(K-1)*log(tau_2_t[u2]*h_pau[[2]][u2])+0.5*sum(log(sigma_alpha[[u2]][v1,]))

        prob2[[v1]][u2] <- expit(w2_cu)

        if (!is.null(s2_cu_oneset[[v1]]) && u2%in%s2_cu_oneset[[v1]]){
          prob2[[v1]][u2] <- 1
        }
      } else{
        mu_alpha[[u2]][v1,] <- (0+(u2==1))*1e10 # basically just to make it c(1,0) for the root.
      }

      # recalculate moments; the ones that matters are the descedant nodes of u1:
      moments_cpp <- get_moments_cpp_eco_alpha_doubletree(
        as.matrix(do.call("rbind",prob2)),
        array(unlist(mu_alpha),c(pL1,K-1,p2)),array(unlist(sigma_alpha),c(pL1,K-1,p2)),
        E_eta,E_eta_sq,ancestors[[2]],leaf_ids_nodes[[2]][[u2]])

      E_eta          <- moments_cpp$E_eta
      E_eta_sq       <- moments_cpp$E_eta_sq
    }
  }

  # update dirich_mat: ---------------------------------------------------
  # the variational parameters for the CSMFs in all the domains:
  for (v2 in 1:pL2){
    for (v1 in 1:pL1){
      dirich_mat[v1,v2] <- dmat[v1,v2] + sum(emat[v_units[[2]]==v2,v1])
    }
  }

  # # intermediate quantities:--------------------------------------------------
  digamma_emat <- as.matrix(sweep(digamma(dirich_mat),MARGIN = 2,digamma(colSums(dirich_mat)),"-"))
  F_array      <- F_doubletree(psi,g_psi,phi,g_phi,X,ind_obs_i,
                               rmat,E_beta,E_beta_sq,E_eta,E_eta_sq,v1_units_NA_replaced,v_units[[2]])
  # update emat:---------------------------------------------------
  if (scenario !="a"){# this means`leaf_ids_units[[1]]$NA_tree1` is not empty:
    emat_update <- update_emat_with_F_doubletree(F_array,rmat,digamma_emat,v_units[[2]])
    emat[leaf_ids_units[[1]]$NA_tree1,] <- as.matrix(emat_update)[leaf_ids_units[[1]]$NA_tree1,]
  }

  # update rmat:-------------------------------------------------------
  # for all subjects their variational probabilities of belonging to each of K classes: step 1b and 1c in Appendix
  if (LD) {rmat <- update_rmat_with_F_doubletree(F_array,emat)}

  # update a_t, b_t: ---------------------------------------------------
  # variational parameters for q_t(rho*_) q_t(rho_cl)
  for (l in 1:Fg1){
    a1_t[l] <- a1[l] + sum(prob1[levels[[1]]==l])
    b1_t[l] <- b1[l] + sum(1-prob1[levels[[1]]==l])
  }

  prob2_mat <- do.call("rbind",prob2)
  for (l in 1:Fg2){
    a2_t[,l] <- a2[,l,drop=FALSE] + rowSums(prob2_mat[,levels[[2]]==l,drop=FALSE])
    b2_t[,l] <- b2[,l,drop=FALSE] + rowSums(1-prob2_mat[,levels[[2]]==l,drop=FALSE])
  }

  # return results:
  make_list(a1_t,b1_t,a2_t,b2_t,emat,
            rmat,dirich_mat,tau_1_t,tau_2_t,
            sigma_gamma,mu_gamma,sigma_alpha,mu_alpha,
            prob1,prob2,
            E_beta_sq,E_eta_sq,E_beta,E_eta)# intermediate calculations useful when calculating ELBO*.
}

#' update hyperparameters
#'
#' NB: argument redundancy may exist
#'
#' @inheritParams update_vi_params_doubletree
#' @param update_hyper Logical, `TRUE` or `FALSE` to indicate
#' whether to update `tau_1` and `tau_2`. This is computed at every iteration
#' in [fit_nlcm_doubletree()]
#' @param E_beta_sq,E_eta_sq,E_beta,E_eta moments computed by [update_vi_params_doubletree()]
#' @param tau_update_levels a numeric vector, specifies which levels of hyperparameters to update
#' @param quiet default to `FALSE`, which prints intermediate updates of hyperparameters
#' @param do_tree1_update `TRUE` to update `mu_gamma`,`sigma_gamma`, `prob1`,`tau_1`
#'
#' @importFrom matrixStats logSumExp
#'
#' @return a list of updated hyperparameters: tau_1,tau_2,psi,g_psi,phi,g_phi,
#'  along with a new ELBO value.
#' @family Internal VI functions
update_hyperparams_doubletree <- function(
  h_pau,
  levels,
  v_units,
  X,n,J,p1,p2,pL1,pL2,Fg1,Fg2, # redundant but streamlined in lcm_tree.
  ind_obs_i,# data and design
  prob1,prob2,mu_gamma,mu_alpha,rmat,emat,dirich_mat,
  sigma_gamma,sigma_alpha,
  tau_1_t,tau_2_t,
  a1_t,b1_t,a2_t,b2_t,
  E_beta_sq,E_eta_sq,E_beta,E_eta,# vi parameters; !!! although this can be calculated from mu's sigma's
  psi,g_psi,phi,g_phi,
  tau_1,tau_2, # hyper-parameters to be updated on a less frequent schedule.
  a1,b1,a2,b2,dmat,K,LD,tau_update_levels,#fixed hyper-parameters not to be update.
  update_hyper,quiet, # called in 'fit_nlcm_tree' to update tau_1 and tau_2 or not.
  do_tree1_update
){
  # calculate some useful quantities:
  v1_units_NA_replaced <- v_units[[1]]
  v1_units_NA_replaced[is.na(v_units[[1]])] <- pL1+1 # this is replacing NA; convenient in C functions.

  prob2_mat <- do.call("rbind",prob2)

  # update local vi parameters:------------------------------------------------
  # the permutation is to match the dimensions in cpp functions
  psi   <- aperm(sqrt(E_beta_sq),c(3,1,2))
  phi   <- aperm(sqrt(E_eta_sq),c(1,3,2))

  g_psi <- g_fun.vec(psi)
  g_phi <- g_fun.vec(phi)

  # update hyper-parameters: tau_1, tau_2:-------------------------------------
  # marginal variational posterior expectation of squared (alpha or gamma):
  expected_ss_alpha <- numeric(p2)
  expected_ss_gamma <- numeric(p1)
  h_pau1 <- h_pau[[1]]
  h_pau2 <- h_pau[[2]]

  for (u in 1:p1){
    expected_ss_gamma[u] <- 1/h_pau1[u]*sum(
      prob1[u]*(sigma_gamma[[u]]+mu_gamma[[u]]^2)+ # J by K
        (1-prob1[u])*matrix(1,nrow=J,ncol=K)*tau_1_t[u]*h_pau1[u])
  }

  for (u in 1:p2){
    expected_ss_alpha[u] <- 1/h_pau2[u]*sum(
      (sigma_alpha[[u]]+(mu_alpha[[u]])^2)*prob2_mat[,u]+# pL1 by K-1; multiplied by pL1
        (1-prob2_mat[,u])*tau_2_t[u]*h_pau2[u]) #pL1 # the sum is pL1*(K-1) elements.
  }

  # update hyperparameters (tau_1, and tau_2) if scheduled: -------------------------------
  # prior is changed if updated.
  if (update_hyper){
    if (do_tree1_update){
      for (l in 1:Fg1){
        if (l %in% tau_update_levels[[1]]){
          tau_1[l]  <- sum(expected_ss_gamma[levels[[1]]==l])/(J*K*sum(levels[[1]]==l))
          cat("> Updated  tau_1; level ",l,":",tau_1[l],". \n")
        }
      }
    }

    for (l in 1:Fg2){
      if (l %in% tau_update_levels[[2]]){
        tau_2[l]  <- sum(expected_ss_alpha[levels[[2]]==l])/(pL1*(K-1)*sum(levels[[2]]==l))
        cat("> Updated  tau_2; level ",l,":",tau_2[l],". \n")
      }
    }
  }

  # update ELBO:
  expected_l_rho2    <- digamma(a2_t) - digamma(a2_t+b2_t)
  expected_l_1m_rho2 <- digamma(b2_t) - digamma(a2_t+b2_t)

  expected_l_rho1    <- digamma(a1_t) - digamma(a1_t+b1_t)
  expected_l_1m_rho1 <- digamma(b1_t) - digamma(a1_t+b1_t)

  # recalculate intermediate quantities because, psi, g_psi, phi, g_phi are now updated.
  digamma_emat <- as.matrix(sweep(digamma(dirich_mat),MARGIN = 2,digamma(colSums(dirich_mat)),"-"))
  F_array      <- F_doubletree(psi,g_psi,phi,g_phi,X,ind_obs_i,
                               rmat,E_beta,E_beta_sq,E_eta,E_eta_sq,v1_units_NA_replaced,v_units[[2]])

  # part 1: E_q(lower bound of joint distribution (all data and unknowns)):
  res1_2_15  <- get_line1_2_15_doubletree(F_array,digamma_emat,rmat,emat,v1_units_NA_replaced,v_units[[2]])
  line_1     <- res1_2_15$res1
  line_2     <- res1_2_15$res2

  line_3     <- - sum(expected_ss_alpha/2/tau_2[levels[[2]]])-sum((K-1)*pL1*log(2*pi*tau_2[levels[[2]]]*h_pau2))/2
  line_4     <- - sum(expected_ss_gamma/2/tau_1[levels[[1]]])-sum(J*K*log(2*pi*tau_1[levels[[1]]]*h_pau1))/2 # this is prior, use hyperparams.

  line_5     <- sum(expected_l_rho2[,levels[[2]],drop=FALSE]*prob2_mat+expected_l_1m_rho2[,levels[[2]],drop=FALSE]*(1-prob2_mat))
  line_6     <- sum(expected_l_rho1[levels[[1]]]*prob1+expected_l_1m_rho1[levels[[1]]]*(1-prob1))

  line_7     <- sum((a2-1)*expected_l_rho2+(b2-1)*expected_l_1m_rho2-mapply(lbeta,c(a2),c(b2)))# here c(a2) c(b2) are to convert a matrix to a vector.
  line_8     <- sum((a1-1)*expected_l_rho1+(b1-1)*expected_l_1m_rho1-mapply(lbeta,a1,b1)) + sum((dmat-1)*digamma_emat) -
    sum(sweep(lgamma(dmat),MARGIN=2,lgamma(colSums(dmat)),"-"))

  # part 2: -E_{q_t}(log(q_t)):
  line_9_10  <- - sum((dirich_mat-1)*digamma_emat) + sum(lgamma(dirich_mat))-sum(lgamma(colSums(dirich_mat)))
  line_11    <- ((K-1)*sum(prob2_mat)*(1+log(2*pi))+ sum(sapply(1:p2,function(u) sum(log(sigma_alpha[[u]])*prob2_mat[,u]))) )/2
  line_12    <-  -((K-1) / 2) * sum(1 - prob2_mat) + ((K-1) / 2) * sum(t(1 - prob2_mat)*log(2 * pi * tau_2_t*h_pau2)) # use tau_2_t because this is VI parameter, happens to be equal to tau_2 numerically.

  line_13    <- (J*K*sum(prob1)*(1+log(2*pi))+sum(sapply(1:p1,function(u) sum(prob1[u]*log(sigma_gamma[[u]])))))/2
  line_14    <- J*K*sum(1-prob1)/2 +J*K*sum(log(2*pi*tau_1_t*h_pau1)*(1-prob1))/2

  line_15    <- res1_2_15$res3

  line_16    <- -1 * (sum(prob2_mat[prob2_mat != 0] * log(prob2_mat[prob2_mat != 0])) +sum((1 - prob2_mat[prob2_mat != 1]) * log(1 - prob2_mat[prob2_mat != 1])))-
    (sum(prob1[prob1 != 0] * log(prob1[prob1 != 0])) +sum((1 - prob1[prob1 != 1]) * log(1 - prob1[prob1 != 1])))

  line_17    <- -1*sum(expected_l_rho2*(a2_t - 1)  +  expected_l_1m_rho2*(b2_t - 1) - mapply(lbeta, c(a2_t), c(b2_t)))# a2_t, b2_t is of pL1 by Fg2 dimension.
  line_18    <- -1*sum(expected_l_rho1*(a1_t - 1)  +  expected_l_1m_rho1*(b1_t - 1) - mapply(lbeta, c(a1_t), c(b1_t)))# a1_t, b1_t is of Fg1 dimension.

  line_vec <- c(line_1, line_2, line_3,line_4,line_5,line_6,line_7,line_8,
                line_9_10, line_11, line_12, line_13,line_14,line_15,line_16,line_17,line_18)

  names(line_vec) <- c("line_1", "line_2", "line_3","line_4","line_5","line_6","line_7","line_8",
                       "line_9_10", "line_11", "line_12", "line_13","line_14","line_15","line_16","line_17","line_18")
  ELBO <- sum(line_vec)

  # return results:
  make_list(ELBO,psi,g_psi,phi,g_phi,tau_1,tau_2)
}
