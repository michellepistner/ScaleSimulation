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
    mutate(low = log2FoldChange - 1.96*lfcSE, 
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
  aldex.fit <- aldex(countdata, as.character(dat$Condition), denom=denom, mc.samples = 1000)
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

###Updated ALDEX clr function

aldex.clr.function <- function( reads, conds, mc.samples=128, denom="all", verbose=FALSE, useMC=FALSE, summarizedExperiment=FALSE, gamma = NULL) {
 ###This is identical to the ALDEx2 package CLR function except lines 223-251
  
  # INPUT
  # The 'reads' data.frame MUST have row
  # and column names that are unique, and
  # looks like the following:
  #
  #              T1a T1b  T2  T3  N1  N2
  #   Gene_00001   0   0   2   0   0   1
  #   Gene_00002  20   8  12   5  19  26
  #   Gene_00003   3   0   2   0   0   0
  #       ... many more rows ...
  #
  # ---------------------------------------------------------------------
  
  # OUTPUT
  # The output returned is a list (x) that contains Monte-Carlo instances of
  # the centre log-ratio transformed values for each sample
  # Access to values
  # sample IDs: names(x)
  # number of features (genes, OTUs): length(x[[1]][,1])
  # number of Monte-Carlo Dirichlet instances: length(x[[1]][1,])
  # feature names: rownames(x[[1]])
  
  # Fully validate and coerce the data into required formats
  # coerce SummarizedExperiment reads into data.frame
  if(summarizedExperiment){
    reads <- data.frame(as.list(assays(reads,withDimnames=TRUE)))
    if (verbose) {
      message("converted SummarizedExperiment read count object into data frame")
    }
  }
  # make sure the conditions vector or matrix is reasonable
  if(missing(conds)){
    
    message("no conditions provided: forcing denom = 'all'")
    message("no conditions provided: forcing conds = 'NA'")
    denom <- "all"
    conds <- rep("NA", ncol(reads))
    
  }
  
  # if a model matrix is supplied, then aldex.effect is not valid
  # force the use of either all for the denominator
  # or
  # the use of a user-supplied denominator
  if(is(conds, "matrix")){
    message("checking for condition length disabled!")
    if(is.vector(denom, mode="integer")){
      message("user-defined denominator used")
    } else if (denom == "all"){
      message("using all features for denominator")
    } else {
      stop("please supply a vector of indices for the denominator")
    }
  }
  
  if(ncol(reads) != length(conds) & !is(conds, "matrix")){
    print(length(conds))
    print(ncol(reads))
    stop("mismatch between number of samples and condition vector")
  }
  
  # make sure that the multicore package is in scope and return if available
  has.BiocParallel <- FALSE
  if ("BiocParallel" %in% rownames(installed.packages()) & useMC){
    message("multicore environment is is OK -- using the BiocParallel package")
    #require(BiocParallel)
    has.BiocParallel <- TRUE
  }
  else {
    message("operating in serial mode")
  }
  
  # make sure that mc.samples is an integer, despite it being a numeric type value
  mc.samples <- as.numeric(as.integer(mc.samples))
  
  #  remove all rows with reads less than the minimum set by minsum
  minsum <- 0
  
  # remove any row in which the sum of the row is 0
  z <- as.numeric(apply(reads, 1, sum))
  reads <- as.data.frame( reads[(which(z > minsum)),]  )
  
  if (verbose) message("removed rows with sums equal to zero")
  
  
  #  SANITY CHECKS ON THE DATA INPUT
  if ( any( round(reads) != reads ) ) stop("not all reads are integers")
  if ( any( reads < 0 ) )             stop("one or more reads are negative")
  
  for ( col in names(reads) ) {
    if ( any( ! is.finite( reads[[col]] ) ) )  stop("one or more reads are not finite")
  }
  
  if ( length(rownames(reads)) == 0 ) stop("rownames(reads) cannot be empty")
  if ( length(colnames(reads)) == 0 ) stop("colnames(reads) cannot be empty")
  
  if ( length(rownames(reads)) != length(unique(rownames(reads))) ) stop ("row names are not unique")
  if ( length(colnames(reads)) != length(unique(colnames(reads))) ) stop ("col names are not unique")
  if ( mc.samples < 128 ) warning("values are unreliable when estimated with so few MC smps")
  
  # add a prior expection to all remaining reads that are 0
  # this should be by a Count Zero Multiplicative approach, but in practice
  # this is not necessary because of the large number of features
  prior <- 0.5
  
  # This extracts the set of features to be used in the geometric mean computation
  # returns a list of features
  if(is.null(gamma)){
    feature.subset <- aldex.set.mode(reads, conds, denom)
    if ( length(feature.subset[[1]]) == 0 ) stop("No low variance, high abundance features in common between conditions\nPlease choose another denomiator.")
  } else{
    feature.subset <- vector()
  }
  
  reads <- reads + prior
  
  if (verbose == TRUE) message("data format is OK")
  
  # ---------------------------------------------------------------------
  # Generate a Monte Carlo instance of the frequencies of each sample via the Dirichlet distribution,
  # returns frequencies for each feature in each sample that are consistent with the
  # feature count observed as a proportion of the total counts per sample given
  # technical variation (i.e. proportions consistent with error observed when resequencing the same library)
  
  nr <- nrow( reads )
  rn <- rownames( reads )
  
  #this returns a list of proportions that are consistent with the number of reads per feature and the
  #total number of reads per sample
  
  # environment test, runs in multicore if possible
  if (has.BiocParallel){
    p <- bplapply( reads ,
                   function(col) {
                     q <- t( rdirichlet( mc.samples, col ) ) ;
                     rownames(q) <- rn ;
                     q })
    names(p) <- names(reads)
  }
  else{
    p <- lapply( reads ,
                 function(col) {
                   q <- t( rdirichlet( mc.samples, col ) ) ;
                   rownames(q) <- rn ; q } )
  }
  
  # sanity check on the data, should never fail
  for ( i in 1:length(p) ) {
    if ( any( ! is.finite( p[[i]] ) ) ) stop("non-finite frequencies estimated")
  }
  
  if (verbose == TRUE) message("dirichlet samples complete")
  
  # ---------------------------------------------------------------------
  # Add scale samples (if desired)
  # Checking the size of the scale samples
  
  if(!is.null(gamma)){
    message("aldex.scaleSim: adjusting samples to reflect scale uncertainty.")
    l2p <- list()
    if(length(gamma) == 1){ ##Add uncertainty around the scale samples
      lambda <- gamma
      scale.samples <-matrix(ncol = mc.samples)
      for(i in 1:length(p)){
        gm_sample <- log(apply(p[[i]],2,gm))
        scale_for_sample <- sapply(gm_sample, FUN = function(mu){stats::rnorm(1, mu, lambda)})
        l2p[[i]] <- sweep(log2(p[[i]]), 2,  log2(exp(scale_for_sample)), "-")
        scale.samples = rbind(scale.samples, scale_for_sample)
      }
      scale.samples <- scale.samples[-1,]
    } else if(length(gamma) >1 & is.null(dim(gamma))){ ##Vector case/scale sim + senstitivity
      stop("A vector was supplied for scale.samples. Please supply only one element.")
    } else{ ##User input of scale samples
      Q <- nrow(gamma)
      N <- ncol(gamma)
      if(Q != ncol(reads) | N != mc.samples){
        stop("Scale samples are of incorrect size!")
      }
      for(i in 1:length(p)){
        l2p[[i]] <- sweep(log2(p[[i]]), 2,  log2(gamma[i,]), "+")
      }
      scale_samples <- gamma
      
    }
    names(l2p) <- names(p)
  }
  
  # ---------------------------------------------------------------------
  # Take the log2 of the frequency and subtract the geometric mean log2 frequency per sample
  # i.e., do a centered logratio transformation as per Aitchison
  
  # apply the function over elements in a list, that contains an array
  if(is.null(gamma)){
    scale_samples <- NULL
    # DEFAULT
    if (is.list(feature.subset)) {
      # ZERO only
      feat.result <- vector("list", length(unique(conds))) # Feature Gmeans
      condition.list <- vector("list", length(unique(conds)))    # list to store conditions
      
      for (i in 1:length(unique(conds)))
      {
        condition.list[[i]] <- which(conds == unique(conds)[i]) # Condition list
        feat.result[[i]] <- lapply( p[condition.list[[i]]], function(m) {
          apply(log2(m), 2, function(x){mean(x[feature.subset[[i]]])})
        })
      }
      set.rev <- unlist(feat.result, recursive=FALSE) # Unlist once to aggregate samples
      p.copy <- p
      for (i in 1:length(set.rev))
      {
        p.copy[[i]] <- as.data.frame(p.copy[[i]])
        p[[i]] <- apply(log2(p.copy[[i]]),1, function(x){ x - (set.rev[[i]])})
        p[[i]] <- t(p[[i]])
      }
      l2p <- p    # Save the set in order to generate the aldex.clr variable
    } else if (is.vector(feature.subset)){
      # Default ALDEx2, iqlr, user defined, lvha
      # denom[1] is put in explicitly for the user-defined denominator case
      if (has.BiocParallel){
        if (denom[1] != "median"){
          l2p <- bplapply( p, function(m) {
            apply( log2(m), 2, function(col) { col - mean(col[feature.subset]) } )
          })
        } else if (denom[1] == "median"){
          l2p <- bplapply( p, function(m) {
            apply( log2(m), 2, function(col) { col - median(col[feature.subset]) } )
          })
        }
        names(l2p) <- names(p)
      }
      else{
        if (denom[1] != "median"){
          l2p <- lapply( p, function(m) {
            apply( log2(m), 2, function(col) { col - mean(col[feature.subset]) } )
          })
        } else if (denom[1] == "median"){
          l2p <- lapply( p, function(m) {
            apply( log2(m), 2, function(col) { col - median(col[feature.subset]) } )
          })
        }
      }
    }  else {
      message("the denominator is not recognized, use a different denominator")
    }
    
    # sanity check on data
    for ( i in 1:length(l2p) ) {
      if ( any( ! is.finite( l2p[[i]] ) ) ) stop("non-finite log-frequencies were unexpectedly computed")
    }
    if (verbose == TRUE) message("transformation complete")
  }
  
  return(new("aldex.clr",reads=reads,mc.samples=mc.samples,conds=conds,denom=feature.subset,verbose=verbose,useMC=useMC,dirichletData=p,analysisData=l2p))
}

run_fakeAldex <- function(dat, n_samples = 2000, gamma = NULL, denom = "all", ...){
  coldata <- dat[,"Condition",drop=F]
  countdata <- t(dat[,-1,drop=F])
  rownames(coldata) <- colnames(countdata) <- paste0("n", 1:ncol(countdata))
  fit <- aldex.clr.function(as.matrix(countdata), as.character(dat$Condition), mc.samples = n_samples, gamma = gamma)
  
  x.tt <- aldex.ttest(fit)
  x.effect <- aldex.effect(fit)
  z <- data.frame(x.effect, x.tt, check.names=F)
  
  return(z)
}


####Scale sim functions for Multinomial Log Normal Model########################
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
  X = rmatnorm(1,mean, diag(nrow(M.star)), V.star)
  Y= C %*% X + M.star
  
  return(Y)
}
 
ssrv.mln<- function(Y = NULL, X = NULL, covariate, upsilon = NULL, Theta = NULL, 
                               Gamma = NULL, Omega = NULL, Xi = NULL, Theta.t = NULL, total_model = "unif", pars = c("Eta", "Lambda", 
                                                                                                   "Sigma"), sample.totals = NULL, sample = NULL, mean_lnorm = NULL, sd_lnorm = NULL, prob = 0.01,...){
  ###Copied from "pibble" function from fido package
  args <- list(...)
  N <- fido:::try_set_dims(c(ncol(Y), ncol(X), args[["N"]]))
  D <- fido:::try_set_dims(c(nrow(Y), nrow(Theta) + 1, nrow(Xi) + 1, ncol(Xi) + 1, args[["D"]]))
  Q <- fido:::try_set_dims(c(nrow(X), ncol(Theta), nrow(Gamma), ncol(Gamma), 
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
  
  calcGradHess <- fido:::args_null("calcGradHess", args, TRUE)
  b1 <- fido:::args_null("b1", args, 0.9)
  b2 <- fido:::args_null("b2", args, 0.99)
  step_size <- fido:::args_null("step_size", args, 0.003)
  epsilon <- fido:::args_null("epsilon", args, 1e-06)
  eps_f <- fido:::args_null("eps_f", args, 1e-10)
  eps_g <- fido:::args_null("eps_g", args, 1e-04)
  max_iter <- fido:::args_null("max_iter", args, 10000)
  verbose <- fido:::args_null("verbose", args, FALSE)
  verbose_rate <- fido:::args_null("verbose_rate", args, 10)
  decomp_method <- fido:::args_null("decomp_method", args, "cholesky")
  eigvalthresh <- fido:::args_null("eigvalthresh", args, 0)
  jitter <- fido:::args_null("jitter", args, 0)
  multDirichletBoot <- fido:::args_null("multDirichletBoot", args, 
                                 -1)
  optim_method <- fido:::args_null("optim_method", args, "lbfgs")
  useSylv <- fido:::args_null("useSylv", args, TRUE)
  ncores <- fido:::args_null("ncores", args, -1)
  seed <- fido:::args_null("seed", args, sample(1:2^15, 1))
  init <- fido:::args_null("init", args, NULL)
  n_samples <- fido:::args_null("n_samples", args, 2000)
  
  if (is.null(init)) 
    init <- fido:::random_pibble_init(Y)
  
  KInv <- chol2inv(chol(Xi))
  AInv <- chol2inv(chol(diag(N) + t(X) %*% Gamma %*% X))
  
  fitc <- optimPibbleCollapsed(Y, upsilon, Theta %*% X, KInv, 
                               AInv, init, n_samples, calcGradHess, b1, b2, step_size, 
                               epsilon, eps_f, eps_g, max_iter, verbose, verbose_rate, 
                               decomp_method, optim_method, eigvalthresh, jitter, multDirichletBoot, 
                               useSylv, ncores, seed)
  ##End of copy from "pibble"; collapsed sampling done
  ##Samples are ALR coordinates
  ##For every total model besides PIM, need to transform back, multiply by the total samples, the transform back
  if(total_model == "pim"){
      tau = matrix(NA, nrow = N, ncol = n_samples)
      D = nrow(Y)
      
      ##Extracting eta.par. Note that F%*% log(W) (== eta) = ALR transform
      lambda.par = fitc$Samples
      
      ##Finding the par and perp components of the priors
      Theta.trans = GG(Theta.t, nrow(Omega))
      Xi.trans = GG(Omega, nrow(Omega))
      
      for(i in 1:n_samples){
        nu.star = upsilon + D - 1
        M.star = Theta.trans$obj.perp %*% X +  Xi.trans$obj.plus %*% solve(Xi.trans$obj.par) %*% ((lambda.par[,,i] - Theta.trans$obj.par %*% X))
        V.star = (diag(N) + t(X) %*% Gamma %*% X) %*% (diag(N) + solve((diag(N) + t(X) %*% Gamma %*% X)) %*% (t(lambda.par[,,i] - Theta.trans$obj.par %*% X) %*% solve(Xi.trans$obj.par) %*% ((lambda.par[,,i] - Theta.trans$obj.par %*% X))))
        Xi.star = Xi.trans$obj.perp - (Xi.trans$obj.plus) %*% solve(Xi.trans$obj.par) %*% t(Xi.trans$obj.plus)
        tau[,i] = t.sampler(nu.star, M.star, Xi.star, V.star)
      }
      
      ##Now, we need to transform the samples back to lambda:
      lambda = array(NA, dim = c(dim(lambda.par)[1] + 1, dim(lambda.par)[2:3]))
      
      ##Defining the transformations that we need
      F <- cbind(diag(D-1), rep(-1,D-1))
      H <- matrix(1,1,D)
      G <- rbind(F,H)
      G.inv <- solve(G)
      for(i in 1:n_samples){
        lambda[,,i] = G.inv %*% rbind(lambda.par[,,i], tau[,i])
      }
      
    } else if(total_model == "flow"){
      tau = matrix(NA, nrow = N, ncol = n_samples)
      if(N != length(sample)){stop("Flow total model does not have the right number of samples!")}
      mean_samp = rep(NA,N)
      sd_samp = rep(NA,N)
      sample.totals$count = log(sample.totals$count)
      
      ##Collapsing by sample to take the average and standard deviation
      for(j in 1:length(sample)){
        tmp.totals = sample.totals %>%
          filter(sample_id == sample[j]) 
        mean_samp[j] = mean(tmp.totals$count)
        sd_samp[j] = sd(tmp.totals$count)
      }
      for(j in 1:length(sample)){
        tau[j,] = rnorm(n_samples, mean_samp[j], sd_samp[j])
      }

    } else if(total_model == "geoMean"){
      tau = matrix(NA, nrow = N, ncol = n_samples)
      if(is.null(sd_lnorm))
        stop("Please provide the standard deviation with total model 'geoMean'.")
      W.par = alrInv_array(fitc$Samples, coords = 1)
      for(j in 1:N){
        geoMean <- apply(W.par[,j,],2, FUN = function(x){exp(mean(log(x)))})
        tau[j,] = rnorm(n_samples, -log(geoMean), sd_lnorm)
      }     
    } else if(total_model == "logNormal"){
      tau = matrix(NA, nrow = N, ncol = n_samples)
      if(is.null(mean_lnorm) | is.null(sd_lnorm))
        stop("Please provide mean and/or standard deviation with total model 'logNormal'.")
      if(length(mean_lnorm) == 1)
        mean_lnorm = rep(mean_lnorm, N)
      if(length(sd_lnorm) == 1)
        sd_lnorm = rep(sd_lnorm, N)
      for(j in 1:N){
        tau[j,] = rnorm(n_samples, mean_lnorm[j], sd_lnorm[j])
      }     
    }

  if(total_model != "pim"){
    ##Combine the samples back together
    lambda.par = alrInv_array(fitc$Samples, coords = 1)
    ##Multiply by the totals
    lambda = array(NA, dim = dim(lambda.par))
    
    for(i in 1:n_samples){
      lambda[,,i] =sweep(log(lambda.par[,,i]), MARGIN=2, tau[,i], `+`)
    }
  }

  grp1 <- which(X[which(rownames(X) == covariate),] == 0)
  grp2 <- which(X[which(rownames(X) == covariate),] == 1)
  
  target.estimator <- apply(lambda, MARGIN = 3, FUN = function(mat, grp1, grp2){rowMeans(mat[,grp1]) - rowMeans(mat[,grp2])}, grp1 = grp1, grp2 = grp2)
  
  mean <- rowMeans(target.estimator)
  low <- apply(target.estimator, 1, FUN = quantile, probs = prob/2)
  high <- apply(target.estimator, 1, FUN = quantile, probs = 1-prob/2)

  return(data.frame(mean = mean, low = low, high = high, category = rownames(Y)))
}#end of function


run_ssrv_mln <- function(dat,
                     upsilon, Gamma,
                     Theta, Omega, Xi, n_samples = 2000, total_model = "unif", sample.totals = NULL, ...){
  coldata <- dat[,"Condition",drop=F]
  countdata <- t(dat[,-1,drop=F])
  rownames(coldata) <- colnames(countdata) <- paste0("n", 1:ncol(countdata))
  fit <- ssrv.mln(as.matrix(countdata), t(model.matrix(~coldata$Condition)),
                             upsilon, Theta, Gamma, Omega, Xi, n_samples = n_samples, total_model = total_model, sample.totals = sample.totals)
  
  return(fit)
}

sig_tram <- function(s){
  filter(s, sign(low)==sign(high))
}

