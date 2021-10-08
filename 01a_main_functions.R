####DeSeq2 functions#############################################################

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

summary_deseq2 <- function(fit, prob=0.05){
  fit %>% 
    as.data.frame() %>% 
    rownames_to_column("category") %>% 
    dplyr::select(category, log2FoldChange, padj, lfcSE) %>% 
    mutate(low = log2FoldChange +1.96*lfcSE, 
           high = log2FoldChange + 1.96*lfcSE) %>% 
    mutate(mean=log2FoldChange)
}

sig_deseq2 <- function(s, pval=0.05){
  filter(s, padj < pval)
}

####Aldex functions##############################################################
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

aldex.ttest <- function(clr, paired.test=FALSE, hist.plot=FALSE, verbose=FALSE) {

  # Use clr conditions slot instead of input
  conditions <- clr@conds
  if(length(unique(conditions)) != 2){
    stop("Please define the aldex.clr object for a vector of two unique 'conditions'.")
  }

  # get dimensions, names, etc from the input data
  smpl.ids <- getSampleIDs(clr)
  feature.names <- getFeatureNames(clr)
  feature.number <- numFeatures(clr)
  mc.instances <- numMCInstances(clr)

  conditions <- as.factor( conditions )
  levels     <- levels( conditions )


  # generate the comparison sets from the condition levels
  levels <- vector( "list", length( levels ) )
  names( levels ) <- levels( conditions )
  sets <- names(levels)
  setAsBinary <- as.numeric(conditions == sets[1])
  setA <- which(conditions == sets[1])
  setB <- which(conditions == sets[2])

  # set up the t-test result containers
  wi.p.matrix <- as.data.frame(matrix(1, nrow = feature.number, ncol = mc.instances))
  wi.BH.matrix <- wi.p.matrix # duplicate result container
  we.p.matrix <- wi.p.matrix # duplicate result container
  we.BH.matrix <- wi.p.matrix # duplicate result container

  # mc.i is the i-th Monte-Carlo instance
  if(verbose) message("running tests for each MC instance:")
  mc.all <- getMonteCarloInstances(clr)
  for(mc.i in 1:mc.instances){

    if(verbose) numTicks <- progress(mc.i, mc.instances, numTicks)

    # generate a matrix of i-th Monte-Carlo instance, columns are samples, rows are features
    t.input <- sapply(mc.all, function(y){y[, mc.i]})

    wi.p.matrix[, mc.i] <- wilcox.fast(t.input, setAsBinary, paired.test)
    wi.BH.matrix[, mc.i] <- p.adjust(wi.p.matrix[, mc.i], method = "BH")

    we.p.matrix[, mc.i] <- t.fast(t.input, setAsBinary, paired.test)
    we.BH.matrix[, mc.i] <- p.adjust(we.p.matrix[, mc.i], method = "BH")
  }

  if(hist.plot == TRUE){

    par(mfrow=c(2,2))
    hist(we.p.matrix[,1], breaks=99, main="Welch's P values Instance 1")
    hist(wi.p.matrix[,1], breaks=99, main="Wilcoxon P values Instance 1")
    hist(we.BH.matrix[,1], breaks=99, main="Welch's BH values Instance 1")
    hist(wi.BH.matrix[,1], breaks=99, main="Wilcoxon BH values Instance 1")
    par(mfrow=c(1,1))
  }

  # get the Expected values of p, q and lfdr
  we.ep <- rowMeans(we.p.matrix) # rowMeans is faster than apply()!!
  we.eBH <- rowMeans(we.BH.matrix)
  wi.ep <- rowMeans(wi.p.matrix)
  wi.eBH <- rowMeans(wi.BH.matrix)

  z <- data.frame(we.ep, we.eBH, wi.ep, wi.eBH)
  rownames(z) <- getFeatureNames(clr)
  return(z)
}

t.fast <- function(data, group, paired){

  grp1 <- group == unique(group)[1]
  grp2 <- group == unique(group)[2]
  n1 <- sum(grp1)
  n2 <- sum(grp2)

  if(paired){

    # Order pairs for the mt.teststat function
    if(n1 != n2) stop("Cannot pair uneven groups.")
    i.1 <- which(grp1)
    i.2 <- which(grp2)
    paired.order <- unlist(lapply(1:length(i.1), function(i) c(i.1[i], i.2[i])))

    t <- multtest::mt.teststat(data[, paired.order], as.numeric(grp1)[paired.order],
                               test = "pairt", nonpara = "n")
    df <- length(i.1) - 1
    return(pt(abs(t), df = df, lower.tail = FALSE) * 2)

  }else{

    t <- multtest::mt.teststat(data, as.numeric(grp1), test = "t", nonpara = "n")
    s1 <- apply(data[, grp1], 1, sd)
    s2 <- apply(data[, grp2], 1, sd)
    df <- ( (s1^2/n1 + s2^2/n2)^2 )/( (s1^2/n1)^2/(n1-1) + (s2^2/n2)^2/(n2-1) )
    return(pt(abs(t), df = df, lower.tail = FALSE) * 2)
  }
}

# wilcox.fast function replaces wilcox.test
#  * runs much faster
#  * uses exact distribution for ties!
#    * this differs from ?wilcox.test
#  * optional paired test
#    * equivalent to wilcox.test(..., correct = FALSE)
#  * uses multtest
wilcox.fast <- function(data, group, paired){

  if(ncol(data) != length(group)) stop("Use rows for feature data.")
  grp1 <- group == unique(group)[1]
  grp2 <- group == unique(group)[2]
  n1 <- sum(grp1)
  n2 <- sum(grp2)

  # Check for ties in i-th Monte-Carlo instance
  data.t <- t(data)
  if(paired){
    anyTies <- any(apply(data.t[grp1, ] - data.t[grp2, ], 2,
                         function(x) length(unique(x))) != ncol(data) / 2)
  }else{
    anyTies <- any(apply(data.t, 2,
                         function(x) length(unique(x))) != ncol(data))
  }

  # Ties trigger slower, safer wilcox.test function
  if(anyTies){
    return(apply(data.t, 2, function(i){
      wilcox.test(x = i[grp1], y = i[grp2], paired = paired, correct = FALSE)$p.value}))
  }

  if(paired){

    if(n1 != n2) stop("Cannot pair uneven groups.")
    data.diff <- data.t[grp1, ] - data.t[grp2, ]
    V <- apply(data.diff, 2, function(x) sum(rank(abs(x))[x > 0]))
    topscore <- (n1 * (n1+1)) / 2
    V.lower <- ifelse(V > topscore / 2, topscore - V, V)
    if(sum(grp1) < 50){ # as per wilcox test, use exact -- ASSUMES NO TIES!!
      V.p <- psignrank(V.lower, n1) * 2
      return(ifelse(V.p > 1, 1, V.p)) # psignrank returns non-zero for W = mean
    }else{ # Use normal approximation
      V.std <- (topscore/2 - V.lower) / sqrt(n1*(n1 + 1) * (2*n1 + 1) / 24) # wilcox.test uses denom = 24
      return(pnorm(V.std, lower.tail = FALSE) * 2)
    }


  }else{

    W.std <- multtest::mt.teststat(data, as.numeric(grp1), test = "wilcoxon")
    if(sum(grp1) < 50 && sum(grp2) < 50){ # as per wilcox test, use exact -- ASSUMES NO TIES!!
      W.var <- sqrt((n1*n2) * (n1+n2+1) / 12)
      W <- abs(W.std) * W.var + (n1*n2) / 2
      W.p <- pwilcox(W - 1, n1, n2, lower.tail = FALSE) * 2
      return(ifelse(W.p > 1, 1, W.p)) # pwilcox returns non-zero for W = mean
    }else{ # Use normal approximation
      return(pnorm(abs(W.std), lower.tail = FALSE) * 2)
    }
  }
}

aldex.kw <- function(clr, useMC=FALSE, verbose=FALSE){

  # Use clr conditions slot instead of input
  conditions <- clr@conds

  # make sure that the multicore package is in scope and return if available
  is.multicore = FALSE

  if ("BiocParallel" %in% rownames(installed.packages()) & useMC){
    message("multicore environment is OK -- using the BiocParallel package")
    #require(BiocParallel)
    is.multicore = TRUE
  }
  else {
    message("operating in serial mode")
  }

  # get dimensions, names, etc from the input data
  smpl.ids <- getSampleIDs(clr)
  feature.names <- getFeatureNames(clr)
  feature.number <- numFeatures(clr)
  mc.instances <- numMCInstances(clr)

  conditions <- as.factor( conditions )
  levels     <- levels( conditions )


  # generate the comparison sets from the condition levels
  levels <- vector( "list", length( levels ) )
  names( levels ) <- levels( conditions )
  sets <- names(levels)

  # set up the glm results containers
  glm.matrix.p <- as.data.frame(matrix(1, nrow = feature.number, ncol = mc.instances))
  glm.matrix.pBH <- glm.matrix.p # duplicate result container
  kw.p.matrix <- glm.matrix.p # duplicate result container
  kw.pBH.matrix <- glm.matrix.p # duplicate result container

  # mc.i is the i-th Monte-Carlo instance
  if(verbose) message("running tests for each MC instance:")
  mc.all <- getMonteCarloInstances(clr)
  for(mc.i in 1:mc.instances){

    if(verbose) numTicks <- progress(mc.i, mc.instances, numTicks)

    # generate a matrix of i-th Monte-Carlo instance, columns are samples, rows are features
    t.input <- sapply(mc.all, function(y){y[, mc.i]})

    # do glms on each feature and make a list of glm outputs
    x <-
      apply(t.input, 1, function(yy) {
        glm(as.numeric(yy) ~ factor(conditions))
      })

    # calculate p-values for generalized linear model
    if(is.multicore == TRUE){
      pps <- bplapply(x, drop1, test = "Chis")
    }else{
      pps <- lapply(x, drop1, test = "Chis")
    }
    glm.matrix.p[, mc.i] <- sapply(pps, function(x){x[[5]][2]})
    glm.matrix.pBH[, mc.i] <- as.numeric(p.adjust(glm.matrix.p[, mc.i], method = "BH"))

    # calculate p-values for Kruskal Wallis test
    kw.p.matrix[, mc.i] <-
      apply(t.input, 1, function(yy){
        kruskal.test(yy, g = factor(conditions))[[3]]
      })
    kw.pBH.matrix[, mc.i] <- as.numeric(p.adjust(kw.p.matrix[, mc.i], method = "BH"))

  }

  # get the Expected values of p, q and lfdr
  glm.ep <- rowMeans(glm.matrix.p) # rowMeans is faster than apply()!!
  glm.eBH <- rowMeans(glm.matrix.pBH)
  kw.ep <- rowMeans(kw.p.matrix)
  kw.eBH <- rowMeans(kw.pBH.matrix)

  z <- data.frame(kw.ep, kw.eBH, glm.ep, glm.eBH)
  rownames(z) <- getFeatureNames(clr)
  return(z)
}


aldex.effect <- function(clr, verbose=TRUE, include.sample.summary=FALSE, useMC=FALSE, CI=FALSE, glm.conds=NULL){

  # Use clr conditions slot instead of input
  if (is.vector(clr@conds)) {
    conditions <- clr@conds
  } else if (is.factor(clr@conds)) {
    if (length(levels(clr@conds) == 2)) {
      conditions <- clr@conds
    }
  } else if (is.matrix(clr@conds)){
    if(is.null(glm.conds)) stop("please provide a binary condition vector")
    conditions <- glm.conds
  } else {
    stop("please check that the conditions parameter for aldex.clr is correct.")
  }

  is.multicore = FALSE

  if ("BiocParallel" %in% rownames(installed.packages()) & useMC==TRUE){
    message("multicore environment is OK -- using the BiocParallel package")
    #require(BiocParallel)
    is.multicore = TRUE
  }
  else {
    if (verbose == TRUE) message("operating in serial mode")
  }

  nr <- numFeatures(clr) # number of features
  rn <- getFeatureNames(clr) # feature names
  # ---------------------------------------------------------------------

  # sanity check to ensure only two conditons passed to this function
  conditions <- as.factor( conditions )
  levels     <- levels( conditions )


  levels <- vector( "list", length( levels ) )
  names( levels ) <- levels( conditions )

  for ( l in levels( conditions ) ) {
    levels[[l]] <- which( conditions == l )
    if ( length( levels[[l]] ) < 2 ) stop("condition level '",l,"' has less than two replicates")
  }

  # end sanity check
  if (verbose == TRUE) message("sanity check complete")

  # Summarize the relative abundance (rab) win and all groups

  rab <- vector( "list", 3 )
  names(rab) <- c( "all", "win", "spl" )
  rab$win <- list()

  #this is the median value across all monte carlo replicates
  cl2p <- NULL
  for ( m in getMonteCarloInstances(clr) ) cl2p <- cbind( cl2p, m )
  rab$all <- t(apply( cl2p, 1, median ))
  rm(cl2p)
  gc()
  if (verbose == TRUE) message("rab.all  complete")

  #this is the median value across all monte carlo replicates per level
  for ( level in levels(conditions) ) {
    cl2p <- NULL
    for ( i in levels[[level]] ) cl2p <- cbind( cl2p, getMonteCarloReplicate(clr,i) )
    rab$win[[level]] <- t(apply( cl2p, 1, median ))
    rm(cl2p)
    gc()
  }
  if (verbose == TRUE) message("rab.win  complete")

  if (is.multicore == TRUE)  rab$spl <- bplapply( getMonteCarloInstances(clr), function(m) { t(apply( m, 1, median )) } )
  if (is.multicore == FALSE) rab$spl <- lapply( getMonteCarloInstances(clr), function(m) { t(apply( m, 1, median )) } )

  if (verbose == TRUE) message("rab of samples complete")

  # ---------------------------------------------------------------------
  # Compute diffs btw and win groups

  l2d <- vector( "list", 2 )
  names( l2d ) <- c( "btw", "win" )
  l2d$win <- list()

  # abs( win-conditions diff ), btw smps
  #this generates a linear sample of the values rather than an exhaustive sample
  for ( level in levels(conditions) ) {
    concat <- NULL
    for ( l1 in sort( levels[[level]] ) ) {
      concat <- cbind(  getMonteCarloReplicate(clr,l1),concat )

    }

    #if the sample is huge, only sample 10000
    if ( ncol(concat) < 10000 ){
      sampl1 <- t(apply(concat, 1, function(x){sample(x, ncol(concat))}))
      sampl2 <- t(apply(concat, 1, function(x){sample(x, ncol(concat))}))
    } else {
      sampl1 <- t(apply(concat, 1, function(x){sample(x, 10000)}))
      sampl2 <- t(apply(concat, 1, function(x){sample(x, 10000)}))
    }
    l2d$win[[level]] <- cbind( l2d$win[[level]] , abs( sampl1 - sampl2 ) )
    rm(sampl1)
    rm(sampl2)
    gc()
  }
  if (verbose == TRUE) message("within sample difference calculated")
  # Handle the case when the groups have different spl sizes
  # get the minimum number of win spl comparisons
  ncol.wanted <- min( sapply( l2d$win, ncol ) )
  # apply multicore paradigm ML
  if (is.multicore == TRUE) l2d$win  <- bplapply( l2d$win, function(arg) { arg[,1:ncol.wanted] } )
  if (is.multicore == FALSE) l2d$win  <- lapply( l2d$win, function(arg) { arg[,1:ncol.wanted] } )

  # btw condition diff (signed)
  #get the btw condition as a random sample rather than exhaustive search
  concatl1 <- NULL
  concatl2 <- NULL
  for( l1 in levels[[1]] ) concatl1 <- cbind( getMonteCarloReplicate(clr,l1),concatl1 )
  for( l2 in levels[[2]] ) concatl2 <- cbind( getMonteCarloReplicate(clr,l2),concatl2 )

  sample.size <- min(ncol(concatl1), ncol(concatl2))

  if ( sample.size < 10000 ){
    smpl1 <- t(apply(concatl1, 1, function(x){sample(x, sample.size)}))
    smpl2 <- t(apply(concatl2, 1, function(x){sample(x, sample.size)}))
  } else {
    smpl1 <- t(apply(concatl1, 1, function(x){sample(x, 10000)}))
    smpl2 <- t(apply(concatl2, 1, function(x){sample(x, 10000)}))
  }
  l2d$btw <- smpl2 - smpl1

  rm(smpl1)
  rm(smpl2)
  gc()
  if (verbose == TRUE) message("between group difference calculated")

  win.max <- matrix( 0 , nrow=nr , ncol=ncol.wanted )
  l2d$effect <- matrix( 0 , nrow=nr , ncol=ncol(l2d$btw) )
  rownames(l2d$effect) <- rn

  ###the number of elements in l2d$btw and l2d$win may leave a remainder when
  #recycling these random vectors. Warnings are suppressed because this is not an issue
  #for this calculation. In fact, any attempt to get rid of this error would
  #decrease our power as one or both vectors would need to be truncated gg 20/06/2013

  options(warn=-1)

  for ( i in 1:nr ) {
    win.max[i,] <- apply( ( rbind( l2d$win[[1]][i,] , l2d$win[[2]][i,] ) ) , 2 , max )
    l2d$effect[i,] <- l2d$btw[i,] / win.max[i,]
  }

  options(warn=0)

  rownames(win.max)   <- rn
  attr(l2d$win,"max") <- win.max
  rm(win.max)

  # ---------------------------------------------------------------------
  # Summarize diffs

  l2s <- vector( "list", 2 )
  names( l2s ) <- c( "btw", "win" )
  l2s$win <- list()

  l2s$btw <- t(apply( l2d$btw, 1, median ))
  l2s$win  <- t(apply( attr(l2d$win,"max"), 1, median ))
  if (verbose == TRUE) message("group summaries calculated")

  if(CI == FALSE) {
    effect  <- t(apply( l2d$effect, 1, median ))
  } else {
    effectlow <- t(apply( l2d$effect, 1, function(x) quantile(x, probs=0.025, names=FALSE) ))
    effecthigh <- t(apply( l2d$effect, 1, function(x) quantile(x, probs=0.975, names=FALSE) ))
    effect  <- t(apply( l2d$effect, 1, median ))
  }
  overlap <- apply( l2d$effect, 1, function(row) { min( aitchison.mean( c( sum( row < 0 ) , sum( row > 0 ) ) + 0.5 ) ) } )
  if (verbose == TRUE) message("effect size calculated")

  # make and fill in the data table
  # i know this is inefficient, but it works and is not a bottleneck
  if(CI == FALSE) {
    rv <- list(
      rab = rab,
      diff = l2s,
      effect = effect,
      overlap = overlap
    )
  } else {
    rv <- list(
      rab = rab,
      diff = l2s,
      effect = effect,
      effectlow = effectlow,
      effecthigh = effecthigh,
      overlap = overlap
    )
  }

  if (verbose == TRUE) message("summarizing output")

  y.rv <- data.frame(t(rv$rab$all))
  colnames(y.rv) <- c("rab.all")
  for(i in names(rv$rab$win)){
    nm <- paste("rab.win", i, sep=".")
    y.rv[,nm] <- data.frame(t(rv$rab$win[[i]]))
  }
  if (include.sample.summary == TRUE){
    for(i in names(rv$rab$spl)){
      nm <- paste("rab.sample", i, sep=".")
      y.rv[,nm] <- data.frame(t(rv$rab$spl[[i]]))
    }

  }
  for(i in names(rv$diff)){
    nm <- paste("diff", i, sep=".")
    y.rv[,nm] <- data.frame(t(rv$diff[[i]]))
  }
  if(CI == FALSE) {
    y.rv[,"effect"] <- data.frame(t(rv$effect))
    y.rv[,"overlap"] <- data.frame(rv$overlap)
  } else {
    y.rv[,"effect"] <- data.frame(t(rv$effect))
    y.rv[,"effect.low"] <- data.frame(t(rv$effectlow))
    y.rv[,"effect.high"] <- data.frame(t(rv$effecthigh))
    y.rv[,"overlap"] <- data.frame(rv$overlap)
  }
  return(y.rv)

}

aitchison.mean <- function( n, log=FALSE ) {

  # Input is a vector of non-negative integer counts.
  # Output is a probability vector of expected frequencies.
  # If log-frequencies are requested, the uninformative subspace is removed.

  n <- round( as.vector( n, mode="numeric" ) )
  if ( any( n < 0 ) ) stop("counts cannot be negative")

  a <- n + 0.5
  sa <- sum(a)

  log.p <- digamma(a) - digamma(sa)
  log.p <- log.p - mean(log.p)

  if ( log ) return(log.p)

  p <- exp( log.p - max(log.p) )
  p <- p / sum(p)
  return(p)
}


####Augmented aldex functions####################################################

fake.aldex <- function(Y = NULL, X = NULL, upsilon = NULL, Theta = NULL, 
                       Gamma = NULL, Omega = NULL, Xi = NULL, total_model = "unif", alpha_total = 1, pars = c("Eta", "Lambda", 
                                                                                                          "Sigma"), test = "t", sample.totals = NULL, n_samples = 2000, mean_lnorm = NULL, sd_lnorm = NULL, ...){
  require(ALDEx2) 
  args <- list(...)
  N <- ncol(Y)
  D <-  nrow(Y)
  Q <- nrow(X)
  if (any(c(N, D, Q) <= 0)) 
    stop("N, D, and Q must all be greater than 0 (D must be greater than 1)")

  ncores <- args_null("ncores", args, -1)
  seed <- args_null("seed", args, sample(1:2^15, 1))
  init <- args_null("init", args, NULL)

  if (is.null(init)) 
    init <- random_pibble_init(Y)
    
  
  
  fitc <- aldex.clr(Y,X,mc.samples = n_samples)
  rownames.hold = row.names(getMonteCarloInstances(fitc)[[1]])
  colnames.hold = colnames(getMonteCarloInstances(fitc)[[1]])
  lambda = getMonteCarloInstances(fitc)
  for(i in 1:length(lambda)){
    lambda[[i]] = clrInv_array(lambda[[i]], 1)
  }
  lambda = array(as.numeric(unlist(lambda)), dim=c(D, n_samples, N))
  
  lambda = aperm(lambda, c(1,3,2))
  
  ##Samples are ALR coordinates
  ##Need to transform back, multiply by the total samples, the transform back
  if(total_model == "unif"){
    if(is.null(sample.totals)){
      print(n_samples)
      tau <- tot_uls(N=N, alpha=alpha_total, n_samples = n_samples)
      print(dim(tau$tau))
    } else{
      tau <- tot_uls(sample=1:length(sample.totals), alpha=alpha_total, w = sample.totals, n_samples = n_samples)
    }}
  if(total_model == "gm"){
    Y.closed <- driver::miniclo_array(Y, parts = 1)
    tau <- matrix(rep(1/apply(Y.closed,2,FUN = gm), n_samples), ncol = n_samples)
    tau =list(tau = tau)
  }
  if(total_model == "logNormal"){
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
    #tau = t(driver::miniclo(t(tau)))
    tau = list(tau = tau)
  }
  ##Transform samples back
  lambda.par = lambda
  ##Multiply by the totals
  lambda = array(NA, dim = dim(lambda.par))
  
  for(i in 1:n_samples){
    lambda[,,i] =sweep(lambda.par[,,i], MARGIN=2, tau$tau[,i]/tau$tau[1,i], `*`)
  }
  if(test == "bayes_glm"){
    collapsed.samples = log(lambda)
    
    
    print(dim(collapsed.samples))
    seed <- seed + sample(1:2^15, 1)
    
    
    ret_mean <- args_null("ret_mean", args, FALSE)
    
    ##Transforming Theta, Gamma, and Xi back to the 
    Xi.t = Omega
    Theta.t = rbind(Theta, rep(0, ncol(Theta)))
    fitu <- uncollapsePibble(collapsed.samples, X, Theta.t, Gamma, 
                             Xi.t, upsilon, ret_mean = ret_mean, ncores = ncores, seed = seed)
    
    out <- list()
    if ("Eta" %in% pars) {
      out[["Eta"]] <- lambda
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
    out$iter <- NA
    out$names_categories <- rownames(Y)
    out$names_samples <- colnames(Y)
    out$names_covariates <- rownames(X)
    out$coord_system <- "base"
    out$alr_base <- D
    out$summary <- NULL
    attr(out, "class") <- c("pibblefit")
    return(out)
  } else if(test == "t"){
    tmp.list = list()
    for(j in 1:dim(lambda)[2]){
      tmp = lambda[,j,]
      ##tmp = apply(tmp, 2, FUN = function(vec) vec/sum(vec))
      ##tmp.list[[j]] = clr_array(tmp, 1)
      tmp.list[[j]] = log(tmp)
      row.names(tmp.list[[j]]) = rownames.hold
    }
    fitc@analysisData = tmp.list
    x.tt <- aldex.ttest(fitc, paired.test = FALSE, hist.plot = FALSE, 
                        verbose = FALSE)
    x.effect <- aldex.effect(fitc, include.sample.summary=FALSE, verbose=FALSE)
    z <- data.frame(x.effect, x.tt, check.names=F)
    return(z)
  } else if(test == "kw"){
    tmp.list = list()
    for(j in 1:dim(lambda)[2]){
      tmp = lambda[,j,]
      tmp = apply(tmp, 2, FUN = function(vec) vec/sum(vec))
      tmp.list[[j]] = clr_array(tmp, 1)
      row.names(tmp.list[[j]]) = rownames.hold
    }
    fitc@analysisData = tmp.list    
    x.tt <- aldex.kw(fitc)
    x.effect <- aldex.effect(fitc, include.sample.summary=FALSE, verbose=FALSE)
    z <- data.frame(x.effect, x.tt, check.names=F)
    return(z)
    
  } else if(test == "lm"){
    tmp.list = list()
    for(j in 1:dim(lambda)[2]){
      tmp = lambda[,j,]
      tmp = apply(tmp, 2, FUN = function(vec) vec/sum(vec))
      tmp.list[[j]] = clr_array(tmp, 1)
      row.names(tmp.list[[j]]) = rownames.hold
    }
    fitc@analysisData = tmp.list    
    x.tt <- aldex.glm(fitc)
    x.effect <- aldex.effect(fitc, include.sample.summary=FALSE, verbose=FALSE)
    z <- data.frame(x.effect, x.tt, check.names=F)
    return(z)
    
  } else(
    return(message("Specified test not supported!"))
  )
  ##Transform samples back
  
  
}#end of function

run_fakeAldex <- function(dat,
                          upsilon, Gamma,
                          Theta, Omega, Xi, n_samples = 2000, total_model = "unif", alpha_total = 0.01, sample.totals = NULL, test=test, mean_lnorm = NULL, sd_lnorm = NULL, ...){
  coldata <- dat[,"Condition",drop=F]
  countdata <- t(dat[,-1,drop=F])
  rownames(coldata) <- colnames(countdata) <- paste0("n", 1:ncol(countdata))
  if(test == "bayes_glm"){
    fit <- fake.aldex(as.matrix(countdata), t(model.matrix(~coldata$Condition)),
                      upsilon, Theta, Gamma, Omega, Xi, n_samples = n_samples, total_model = total_model, alpha_total = alpha_total, sample.totals = sample.totals, test = test, mean_lnorm = mean_lnorm, sd_lnorm = sd_lnorm)
  } else{
    fit <- fake.aldex(as.matrix(countdata), as.character(dat$Condition),
                      upsilon, Theta, Gamma, Omega, Xi, n_samples = n_samples, total_model = total_model, alpha_total = alpha_total, sample.totals = sample.totals, test = test, mean_lnorm = mean_lnorm, sd_lnorm = sd_lnorm)
  }
  
  
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
  Sigma = rinvwishart(nu.star, Xi.star)
  C = t(chol(Sigma))
  X = rmatnorm(1,M.star,diag(nrow(M.star)), V.star)
  Y= C %*% X
  
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
  if(is.null(total_model)){
    tau <- tau_fit
  } else{
    if(total_model == "unif"){
      if(is.null(w)){
        tau <- tot_uls(N=N, alpha=alpha_total)
      } else{
        tau <- tot_uls(sample=sample, alpha=alpha_total, w = w)
      }}
    if(total_model == "unif_half"){
      tau <- tot_half(N=N, alpha=alpha_total)
      tau = list(tau=tau)
    }
    if(total_model == "gm"){
      Y.closed <- driver::miniclo_array(Y, parts = 1)
      tau <- matrix(rep(1/apply(Y.closed,2,FUN = gm), n_samples), ncol = n_samples)
      tau =list(tau = tau)
    } 
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
    }
    if(total_model == "flow"){
      tau = matrix(NA, nrow = N, ncol = n_samples)
      if(N != length(sample)){stop("Flow total model does not have the right number of samples!")}
      for(j in 1:length(sample)){
        tmp.totals = sample.totals %>%
          filter(sample_id == sample[j]) 
        mean_samp = mean(tmp.totals$count)
        sd_samp = sd(tmp.totals$count)
        tau[j,] = rnorm(n_samples, mean_samp, sd_samp)
      }
      tau = t(driver::miniclo(t(tau)))
      tau = list(tau = tau)
    }
    if(total_model == "logNormal"){
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
      #tau = t(driver::miniclo(t(tau)))
      tau = list(tau = tau)
    }
    if(total_model == "logNormalpairs"){
      tau = matrix(NA, nrow = N, ncol = n_samples)
      if(is.null(mean_lnorm) | is.null(sd_lnorm))
        stop("Please provide mean and/or standard deviation with total model 'logNormal'.")
      tau[(N/2+1):N,] = 1 ## runif((N-(N/2))*n_samples)
      for(j in 1:(N/2)){
        tau[j,] = exp(rnorm(n_samples, mean_lnorm, sd_lnorm)) * tau[(N/2+j),]
      }
      #tau = t(driver::miniclo(t(tau)))
      tau = list(tau = tau)
    }}

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
    
    
    # FF = cbind(diag(D-1), rep(-1,D-1))
    # HH = matrix(1,1,D)
    # G.mat = rbind(FF,HH)
    # lambda.par_tmp = array(NA, dim = c(dim(lambda.par)[1] - 1, dim(lambda.par)[2:3]))
    # lambda = array(NA, dim = dim(lambda.par)) #Always needed
    # for(i in 1:n_samples){
    #   lambda.par_tmp[,,i] = FF %*% lambda.par[,,i]
    # }
    # for(i in 1:n_samples){
    #   lambda[,,i] = solve(G.mat) %*% rbind(lambda.par_tmp[,,i], tau$tau[,i])
    # }
    # collapsed.samples = lambda
  } else{
    lambda.par = alrInv_array(fitc$Samples, coords = 1)
    ##Multiply by the totals
    lambda = array(NA, dim = dim(lambda.par))
    
    for(i in 1:n_samples){
      lambda[,,i] =sweep(lambda.par[,,i], MARGIN=2, tau$tau[,i]/tau$tau[1,i], `*`)
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
  }
  else {
    ret_mean <- args_null("ret_mean", args, FALSE)
  }
  
  ##Transforming Theta, Gamma, and Xi back to the 
  Xi.t = Omega
  Theta.t = rbind(Theta, rep(0, ncol(Theta)))
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


###Total model
totulsfit = function (N, iter, tau, names_samples, log_scale, K, w = NULL, 
                      sample = NULL, alpha = NULL) 
{
  tf <- new_totfit(N, iter, tau, names_samples, log_scale = log_scale, 
                   K = K, w = w, sample = sample, alpha = alpha)
  return(tf)
}

new_totfit <- function(N, iter, tau, names_samples, ...){
  args <- list(...)
  contents <- list(N=N, iter=iter, tau=tau, names_samples=names_samples)
  contents <- c(contents, args)
  out <- structure(contents,  class=c("totfit"))
  return(out)
}

tot_uls = function (w = NULL, sample = NULL, N = NULL, alpha = NULL, n_samples = 5000, 
                    log_scale = FALSE, use_names = TRUE) 
{
  K <- length(w)
  if (!is.null(w) & is.null(sample)) 
    stop("sample must be non-null if w is passed")
  if (!is.null(w) & !is.null(sample)) {
    if (log_scale) {
      w <- log(w)
    }
    if (!is.factor(sample)) 
      sample <- as.factor(sample)
    if (!is.null(N)) 
      stopifnot(N == length(levels(sample)))
    if (length(sample) != length(w)) 
      stop("sample and w have different lengths")
    N <- length(levels(sample))
    m <- split(w, sample)
    m <- unlist(lapply(m, mean))
    m <- m[levels(sample)]
    if (log_scale) {
      w <- exp(w)
      m <- exp(m)
    }
  }
  else {
    stopifnot(!is.null(N))
    m <- NULL
  }
  tau <- driver::rUnifSphere(n_samples, N - 1, alpha, shell_only = FALSE)
  if (!is.null(w)) {
    m <- ilr(m)
    tau <- sweep(tau, 1, m, FUN = `+`)
  }
  tau <- t(ilrInv(t(tau)))
  character.levels <- is_charachter_vector(sample)
  names_samples <- switch(as.integer(character.levels) + 1, 
                          NULL, levels(sample))
  out <- totulsfit(N = as.integer(N), K = as.integer(K), iter = as.integer(n_samples), 
                   w = w, sample = sample, alpha = alpha, tau = tau, names_samples = names_samples, 
                   log_scale = log_scale)

  return(out)
}

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

tot_half = function (N = NULL, alpha = NULL, alpha.lower = 0, n_samples = 2000, samples = NULL, shell_only = FALSE) 
{
  
  z <- matrix(rnorm(N*n_samples), N, n_samples)
  z[7:11,] <- abs(z[7:11,])
  sf <- sqrt(colSums(z^2))
  if (shell_only) {
    ru <- radius
  }
  else {
    ru <- alpha * runif(n_samples)^(1/N)
  }
  tau <- sweep(z, 2, ru/sf, FUN = "*")
  
  tau <- abs(tau)
  tau <- t(ilrInv(t(tau)))
  return(tau)
}

