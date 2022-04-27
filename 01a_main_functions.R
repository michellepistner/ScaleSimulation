####DeSeq2 functions#############################################################

###DESeq2 master function for simulation study
###Input: OTU table with condition appended to it
###Output: A DESeq2 results function
run_deseq2 <- function(dat){
  coldata <- dat[,"Condition",drop=F]
  countdata <- t(dat[,-1,drop=F])
  rownames(coldata) <- colnames(countdata) <- paste0("n", 1:ncol(countdata))
  dds <- DESeqDataSetFromMatrix(countData=countdata,
                                colData=coldata,
                                design = ~Condition)
  dds <- DESeq(dds)
  res <- results(dds)
  return(res)
}

##Summarizes DESeq2 results function
summary_deseq2 <- function(fit, prob=0.05){
  fit %>% 
    as.data.frame() %>% 
    rownames_to_column("category") %>% 
    dplyr::select(category, log2FoldChange, padj, lfcSE) %>% 
    mutate(low = log2FoldChange +1.96*lfcSE, 
           high = log2FoldChange + 1.96*lfcSE) %>% 
    mutate(mean=log2FoldChange)
}

##Filters DESEq summary object to those taxa with p-values less than a threshold
sig_deseq2 <- function(s, pval=0.05){
  filter(s, padj < pval)
}

####Aldex functions##############################################################
##Master function to run ALDEx2 for simulated data
##Input: OTU table with condition column appended to it
##Output: An aldex object
run_aldex2 <- function(dat, denom="all"){
  countdata <- t(dat[,-1,drop=F])
  colnames(countdata) <- paste0("n", 1:ncol(countdata))
  aldex.fit <- aldex(countdata, as.character(dat$Condition), denom=denom)
  return(aldex.fit)
}

summary_aldex2 <- function(fit){
  fit %>% 
    as.data.frame() %>% 
    rownames_to_column("category") %>% 
    dplyr::select(category, effect, wi.eBH) %>% 
    mutate(padj=wi.eBH) %>% 
    mutate(mean=effect) %>% 
    mutate(low=NA, high=NA)
}

sig_aldex2 <- function(s, pval=0.05){
  filter(s, padj < pval)
}

####Augmented aldex functions####################################################

run_fakeAldex <- function(dat, n_samples = 2000,  test=test, scale.samples = NULL, ...){
  coldata <- dat[,"Condition",drop=F]
  countdata <- t(dat[,-1,drop=F])
  rownames(coldata) <- colnames(countdata) <- paste0("n", 1:ncol(countdata))
  fit <- aldex(as.matrix(countdata), as.character(dat$Condition), mc.samples = n_samples,test = test, gamma = scale.samples)

  return(fit)
}


####Supplementation functions####################################################
GG <- function(obj, D){
  if(is.null(nrow(obj)) | nrow(obj) != ncol(obj)){
    F = cbind(diag(D-1), rep(-1,D-1))
    H = matrix(1,1,D)
    
    obj.par = F%*% obj
    obj.perp = H %*% obj
    return(list(obj.par = obj.par, obj.perp = obj.perp))
  } else{
    F = cbind(diag(D-1), rep(-1,D-1))
    H = matrix(1,1,D)
    
    obj.par = F%*% obj %*% t(F)
    obj.perp = H %*% obj %*% t(H)
    obj.plus = H %*% obj %*% t(F)
    return(list(obj.par = obj.par, obj.perp = obj.perp, obj.plus = obj.plus))
  }
}

t.sampler <- function(nu.star, M.star, Xi.star, V.star){
  Sigma = rinvwishart(nu.star + 2*nrow(Xi.star), Xi.star)
  C = t(chol(Sigma))
  mean = matrix(0, nrow = nrow(M.star), ncol = ncol(M.star))
  X = rmatnorm(1,mean,diag(nrow(M.star)), V.star)
  Y= C %*% X + M.star
  
  return(Y)
}



supplementation.mln<- function(Y = NULL, X = NULL, upsilon = NULL, Theta = NULL, 
                               Gamma = NULL, Omega = NULL, Xi = NULL, total_model = "unif", pars = c("Eta", "Lambda", 
                                                                                                     "Sigma"), sample.totals = NULL, alpha_total = 1, tau_fit = NULL, sample = NULL, w = NULL, mean_lnorm = NULL, sd_lnorm = NULL, ...){
  args <- list(...)
  N <- try_set_dims(c(ncol(Y), ncol(X), args[["N"]]))
  D <- try_set_dims(c(nrow(Y), nrow(Theta) + 1, nrow(Xi) + 1, ncol(Xi) + 1, args[["D"]]))
  Q <- try_set_dims(c(nrow(X), ncol(Theta), nrow(Gamma), ncol(Gamma), 
                      args[["Q"]]))
  if (any(c(N, D, Q) <= 0)) 
    stop("N, D, and Q must all be greater than 0 (D must be greater than 1)")
  if (D <= 1) 
    stop("D must be greater than 1")
  if (is.null(upsilon)) 
    upsilon <- D + 3
  if (is.null(Theta)) 
    Theta <- matrix(0, D - 1, Q)
  if (is.null(Gamma)) 
    Gamma <- diag(Q)
  if (is.null(Xi)) {
    Xi <- matrix(0.5, D - 1, D - 1)
    diag(Xi) <- 1
    Xi <- Xi * (upsilon - D)
  }
  
  calcGradHess <- args_null("calcGradHess", args, TRUE)
  b1 <- args_null("b1", args, 0.9)
  b2 <- args_null("b2", args, 0.99)
  step_size <- args_null("step_size", args, 0.003)
  epsilon <- args_null("epsilon", args, 1e-06)
  eps_f <- args_null("eps_f", args, 1e-10)
  eps_g <- args_null("eps_g", args, 1e-04)
  max_iter <- args_null("max_iter", args, 10000)
  verbose <- args_null("verbose", args, FALSE)
  verbose_rate <- args_null("verbose_rate", args, 10)
  decomp_method <- args_null("decomp_method", args, "cholesky")
  eigvalthresh <- args_null("eigvalthresh", args, 0)
  jitter <- args_null("jitter", args, 0)
  multDirichletBoot <- args_null("multDirichletBoot", args, 
                                 -1)
  optim_method <- args_null("optim_method", args, "lbfgs")
  useSylv <- args_null("useSylv", args, TRUE)
  ncores <- args_null("ncores", args, -1)
  seed <- args_null("seed", args, sample(1:2^15, 1))
  init <- args_null("init", args, NULL)
  n_samples <- args_null("n_samples", args, 2000)
  
  if (is.null(init)) 
    init <- random_pibble_init(Y)
  
  KInv <- chol2inv(chol(Xi))
  AInv <- chol2inv(chol(diag(N) + t(X) %*% Gamma %*% X))
  
  fitc <- optimPibbleCollapsed(Y, upsilon, Theta %*% X, KInv, 
                               AInv, init, n_samples, calcGradHess, b1, b2, step_size, 
                               epsilon, eps_f, eps_g, max_iter, verbose, verbose_rate, 
                               decomp_method, optim_method, eigvalthresh, jitter, multDirichletBoot, 
                               useSylv, ncores, seed)
  
  ##Samples are ALR coordinates
  ##Need to transform back, multiply by the total samples, the transform back
  if(total_model == "pim"){
      tau = matrix(NA, nrow = N, ncol = n_samples)
      D = nrow(Y)
      
      lambda.par = alrInv_array(fitc$Samples, coords = 1)
      FF = cbind(diag(D-1), rep(-1,D-1))
      lambda.par_tmp = array(NA, dim = c(dim(lambda.par)[1] - 1, dim(lambda.par)[2:3]))
      for(i in 1:n_samples){
        lambda.par_tmp[,,i] = FF %*% lambda.par[,,i]
      }
      
      Xi.t = Omega
      Theta.t = rbind(Theta, rep(0, ncol(Theta)))
      
      Theta.trans = GG(Theta.t, nrow(Xi.t))
      Xi.trans = GG(Xi.t, nrow(Xi.t))
      
      FF = cbind(diag(D-1), rep(-1,D-1))
      for(i in 1:n_samples){
        nu.star = upsilon + D
        M.star = Theta.trans$obj.perp %*% X +  Xi.trans$obj.plus %*% solve(Xi.trans$obj.par) %*% ((lambda.par_tmp[,,i] - Theta.trans$obj.par %*% X))
        V.star = (diag(N) + t(X) %*% Gamma %*% X) %*% (diag(N) + solve((diag(N) + t(X) %*% Gamma %*% X)) %*% (t(lambda.par_tmp[,,i] - Theta.trans$obj.par %*% X) %*% solve(Xi.trans$obj.par) %*% ((lambda.par_tmp[,,i] - Theta.trans$obj.par %*% X))))
        Xi.star = Xi.trans$obj.perp - (Xi.trans$obj.plus) %*% solve(Xi.trans$obj.par) %*% t(Xi.trans$obj.plus)
        tau[,i] = t.sampler(nu.star, M.star, Xi.star, V.star)
      }
      tau = list(tau = tau)
    } else if(total_model == "flow"){
      tau = matrix(NA, nrow = N, ncol = n_samples)
      if(N != length(sample)){stop("Flow total model does not have the right number of samples!")}
      mean_samp = rep(NA,N)
      sd_samp = rep(NA,N)
      sd_log = rep(NA,N)
      sample.totals$count = sample.totals$count/1e9
      for(j in 1:length(sample)){
        tmp.totals = sample.totals %>%
          filter(sample_id == sample[j]) 
        mean_samp[j] = mean(tmp.totals$count)
        sd_samp[j] = sd(tmp.totals$count)
        sd_log[j] = sqrt(sd_samp[j]^2/(mean_samp[j]^2))
      }
      for(j in 1:length(sample)){
        tau[j,] = rlnorm(n_samples, log(mean_samp[j]), sd_log[j])
      }

      tau = list(tau = tau)
    } else if(total_model == "logNormal"){
      tau = matrix(NA, nrow = N, ncol = n_samples)
      if(is.null(mean_lnorm) | is.null(sd_lnorm))
        stop("Please provide mean and/or standard deviation with total model 'logNormal'.")
      if(length(mean_lnorm) == 1)
        mean_lnorm = rep(mean_lnorm, N)
      if(length(sd_lnorm) == 1)
        sd_lnorm = rep(sd_lnorm, N)
      for(j in 1:N){
        tau[j,] = rlnorm(n_samples, mean_lnorm[j], sd_lnorm[j])
      }     
      tau = list(tau = tau)
    }

  ##Transform samples back
  if(total_model == "pim"){
    W.par = alrInv_array(fitc$Samples, coords = 1)
    lambda.par = log(W.par)
    lambda = array(NA, dim = dim(lambda.par)) #Always needed
    for(i in 1:n_samples){
      lambda[,,i] =sweep(lambda.par[,,i], MARGIN=2, tau$tau[,i], `+`)
    }
    ##Transform samples back
    collapsed.samples = lambda
    
  } else{
    lambda.par = alrInv_array(fitc$Samples, coords = 1)
    ##Multiply by the totals
    lambda = array(NA, dim = dim(lambda.par))
    
    for(i in 1:n_samples){
      lambda[,,i] =sweep(lambda.par[,,i], MARGIN=2, tau$tau[,i], `*`)
    }
    ##Transform samples back
    collapsed.samples = log(lambda)
  }
  

  print(dim(collapsed.samples))
  seed <- seed + sample(1:2^15, 1)
  if (is.null(fitc$Samples)) {
    fitc$Samples <- add_array_dim(fitc$Pars, 3)
    ret_mean <- args_null("ret_mean", args, TRUE)
    if (ret_mean && n_samples > 0) {
      warning("Laplace Approximation Failed, using MAP estimate of eta", 
              " to obtain Posterior mean of Lambda and Sigma", 
              " (i.e., not sampling from posterior distribution of Lambda or Sigma)")
    }
    if (!ret_mean && n_samples > 0) {
      warning("Laplace Approximation Failed, using MAP estimate of eta", 
              "but ret_mean was manually specified as FALSE so sampling", 
              "from posterior of Lambda and Sigma rather than using posterior mean")
    }
  } else {
    ret_mean <- args_null("ret_mean", args, FALSE)
  }
  
  ##Transforming Theta, Gamma, and Xi back to the 
  Xi.t = Omega
  Theta.t = alrInv_array(Theta, D, 1)
    rbind(Theta, rep(0, ncol(Theta)))
  fitu <- uncollapsePibble(collapsed.samples, X, Theta.t, Gamma, 
                           Xi.t, upsilon, ret_mean = ret_mean, ncores = ncores, seed = seed)
  
  out <- list()
  if ("Eta" %in% pars) {
    out[["Eta"]] <- collapsed.samples
  }
  if ("Lambda" %in% pars) {
    out[["Lambda"]] <- fitu$Lambda
  }
  if ("Sigma" %in% pars) {
    out[["Sigma"]] <- fitu$Sigma
  }
  out$N <- N
  out$Q <- Q
  out$D <- D
  out$Y <- Y
  out$upsilon <- upsilon
  out$Theta <- Theta
  out$X <- X
  out$Xi <- Xi
  out$Gamma <- Gamma
  out$init <- init
  out$iter <- dim(fitc$Samples)[3]
  out$names_categories <- rownames(Y)
  out$names_samples <- colnames(Y)
  out$names_covariates <- rownames(X)
  out$coord_system <- "base"
  out$alr_base <- D
  out$summary <- NULL
  attr(out, "class") <- c("pibblefit")
  return(out)
  
}#end of function


run_tram <- function(dat,
                     upsilon, Gamma,
                     Theta, Omega, Xi, n_samples = 2000, total_model = "unif", alpha_total = 0.01, sample.totals = NULL, ...){
  coldata <- dat[,"Condition",drop=F]
  countdata <- t(dat[,-1,drop=F])
  rownames(coldata) <- colnames(countdata) <- paste0("n", 1:ncol(countdata))
  fit <- supplementation.mln(as.matrix(countdata), t(model.matrix(~coldata$Condition)),
                             upsilon, Theta, Gamma, Omega, Xi, n_samples = n_samples, total_model = total_model, alpha_total = alpha_total, sample.totals = sample.totals)
  
  return(fit)
}

summary_tram <- function(fit, prob=.9){
  prob1 = .80
  prob2 = .99
  df.results = tibble(.rows = fit$D * fit$Q) %>%
    mutate("category" = rep(fit$names_categories, fit$Q)) %>%
    mutate("covariate" = rep(fit$names_covariates, each = fit$D)) %>%
    mutate("low" = c(apply(fit$Lambda, c(1,2), min))) %>%
    mutate("pLow" = c(apply(fit$Lambda, c(1,2), quantile, prob = (1-prob)/2, na.rm = TRUE))) %>%
    mutate("pLow80" = c(apply(fit$Lambda, c(1,2), quantile, prob = (1-prob1)/2, na.rm = TRUE))) %>%
    mutate("pLow99" = c(apply(fit$Lambda, c(1,2), quantile, prob = (1-prob2)/2, na.rm = TRUE))) %>%
    mutate("DE" = c(apply(fit$Lambda, c(1,2), mean, na.rm = TRUE))) %>%
    mutate("pHigh" = c(apply(fit$Lambda, c(1,2), quantile, prob = prob+(1-prob)/2, na.rm = TRUE))) %>%
    mutate("pHigh80" = c(apply(fit$Lambda, c(1,2), quantile, prob = prob1+(1-prob1)/2, na.rm = TRUE))) %>%
    mutate("pHigh99" = c(apply(fit$Lambda, c(1,2), quantile, prob = prob2+(1-prob2)/2, na.rm = TRUE))) %>%
    mutate("high" = c(apply(fit$Lambda, c(1,2), max, na.rm = TRUE))) %>%
    filter(covariate != "(Intercept)")
  
  return(df.results)
}

sig_tram <- function(s){
  filter(s, sign(pLow)==sign(pHigh))
}

