###Loading libraries
library(plyr)
library(dplyr)
library(phyloseq)
library(fido)
library(driver)
library(tidyverse)
library(DESeq2)
library(ALDEx2)
library(gghighlight)
library(cowplot)
library(ggplot2)
library(magrittr)
library(stringi)
library(directlabels)
library(LaplacesDemon)
library(MASS)
library(matrixNormal)
library(ggpattern)
library(latex2exp)
set.seed(1993)

###Setting the data parameters for all simulations
d <- c(4000, 4000, 4000, 4000, 4000, 400,400,400,400,4000,400,500,500,500,400,400,400,400,400,100,400, # Pre
       4000, 4000, 3000, 2000, 4000, 400,400,400,400,4000,400,500,500,500,200,400,400,400,400,100,100) # Post
dd = length(d)/2
truth1 <- !(d[1:dd]==d[(dd+1):(2*dd)])##testing if the mean is different



aldexDeseq_analysis <- function(d, n, seq.depth, pval = 0.05, prob = .99){
  dd <- length(d)/2
  truth1 <- !(d[1:dd]==d[(dd+1):(2*dd)])##testing if the mean is different
  truth2 <- (1:dd)[truth1]##Locations of the differences
  dat <- create_true_abundances(d, n=n)
  rdat <- resample_data(dat, seq.depth=seq.depth)
  
  dd <- ncol(rdat)-1
  dfit <- run_deseq2(rdat)
  dfit <- summary_deseq2(dfit)
  dfit <- append_sig(dfit, function(x) sig_deseq2(x, pval=pval))
  
  afit <- run_aldex2(rdat)
  afit <- summary_aldex2(afit)
  afit <- append_sig(afit, function(x) sig_aldex2(x, pval=pval))
  
  upsilon = dd + 10
  Gamma = diag(c(5,1))
  Theta = matrix(0, dd-1, 2)
  Omega = diag(dd)
  G = cbind(diag(dd-1), -1)
  Xi = G%*%Omega%*%t(G)
  sfit = run_fakeAldex(rdat,
                       upsilon = upsilon, Gamma = Gamma,
                       Theta = Theta, Omega = Omega, Xi = Xi, n_samples = 5000, total_model = "logNormal", mean_lnorm = 0, sd_lnorm = 1, test="bayes_glm")
  
  sfit <- summary_tram(sfit, prob=prob)
  sfit <- append_sig(sfit, sig_tram)
  
  rrs <- list(dat=rdat, dfit=dfit, afit=afit)
  
  p1 <- plot_count(dat)
  p2 <- plot_sig2(rrs, truth=truth2)
  p1 <- p1+theme(axis.title.x = element_blank(), 
                 axis.text.x = element_blank(), 
                 axis.ticks.x=element_blank(),
                 text = element_text(size=16))
  p <- plot_grid(p1, p2, nrow=2, align="v", rel_heights=c(1.7, 1))
  p
  
  fp.deseq = sum(dfit$sig == TRUE & truth1 == FALSE)
  tp.deseq = sum(dfit$sig == TRUE & truth1 == TRUE)
  
  tp.sfit = sum(sfit$sig == TRUE & truth1 == TRUE)
  fp.sfit = sum(sfit$sig == TRUE & truth1 == FALSE)
  
  fp.aldex = sum(afit$sig == TRUE & truth1 == FALSE)
  tp.aldex = sum(afit$sig == TRUE & truth1 == TRUE)
  
  
  return(list("pval_deseq" = dfit$padj, "fp_deseq" = fp.deseq, "tp_deseq" = tp.deseq, "fp_sfit" = fp.sfit, "tp_sfit" = tp.sfit, "pval_aldex" = afit$padj, "fp_aldex" = fp.aldex, "tp_aldex" = tp.aldex))
  
  #ggsave(p1, file.path("Results","aldex_deseq_failures.pdf"))
}# end of function

model.names <- c("dfit"="DESeq2",
                 "afit"="ALDEx2",
                 "sfit"= "Scale Simulation (Relaxed)")
model.name.levels <- c("DESeq2", "ALDEx2", "Scale Simulation (Relaxed)")


vals = c(5,10,25,50,75,100, 125, 150, 200)

fdr.aldex = rep(NA, length(vals))
fdr.deseq = rep(NA, length(vals))
fdr.sfit = rep(NA, length(vals))

pval.aldex = matrix(NA, ncol = length(d)/2, nrow = length(vals))
pval.deseq = matrix(NA, ncol = length(d)/2, nrow = length(vals))
pval.sfit = matrix(NA,ncol = length(d)/2, nrow = length(vals))


for(i in 1:length(vals)){
  res = aldexDeseq_analysis(d, vals[i], 5000)
  
  fdr.aldex[i] = res$fp_aldex/(res$fp_aldex + res$tp_aldex)
  fdr.deseq[i] = res$fp_deseq/(res$fp_deseq + res$tp_deseq)
  fdr.sfit[i] = res$fp_sfit/(res$fp_sfit + res$tp_sfit)
  fdr.sfit[i] = ifelse((res$fp_sfit + res$tp_sfit) == 0, 0, fdr.sfit[i])
  
  pval.aldex[i,] = res$pval_aldex
  pval.deseq[i,] = res$pval_deseq

  print(i)
}

graph.aldex = data.frame()
for(i in 1:ncol(pval.aldex)){
  tmp = cbind(vals, pval.aldex[,i], rep(paste0("Taxa ", i), length(vals)), rep(truth1[i], length(vals)))
  graph.aldex = rbind(graph.aldex, tmp)
}

names(graph.aldex) = c("vals", "pval", "taxa", "sig")
graph.aldex$pval = as.numeric(graph.aldex$pval)
graph.aldex$vals = as.numeric(graph.aldex$vals)

fdr.all = data.frame(vals = rep(vals,2), fdr = c(c(fdr.aldex),c(fdr.deseq)), method = rep(c("ALDEx2", "DESeq2"), each = length(vals)))
fdr.all$method = as.factor(fdr.all$method)

ggplot(fdr.all, aes(x=vals, y=fdr, color=method, fill = method, linetype = method)) +
  geom_line(alpha = 1, lwd = 1.1) +
  theme_bw()+
  xlim(c(0,200)) +
  xlab("Sample Size") +
  ylab("False Discovery Rate") + 
  scale_color_manual(values=c("#5e3c99", "#e66101")) + 
  scale_linetype_manual(values = c("dotted", "twodash")) +
  theme(text=element_text(size=14)) + 
  geom_hline(yintercept=17/21, linetype="dashed", color = "grey") +
  theme(legend.title = element_blank())
ggsave(file.path("results", "unacknowledged_bias.pdf"), height=3, width=6)

fdr.all = data.frame(vals = rep(vals,3), fdr = c(c(fdr.aldex),c(fdr.deseq), c(fdr.sfit)), method = rep(c("ALDEx2", "DESeq2", "Scale Sim. (Relaxed)"), each = length(vals)))
fdr.all$method = as.factor(fdr.all$method)

ggplot(fdr.all, aes(x=vals, y=fdr, color=method, fill = method, linetype = method)) +
  geom_line(alpha = 1, lwd = 1.1) +
  theme_bw()+
  xlim(c(0,200)) +
  xlab("Sample Size") +
  ylab("False Discovery Rate") + 
  scale_color_manual(values=c("#AF4BCE", "#EA7369", "#2b83ba")) + 
  scale_linetype_manual(values = c("dotted", "twodash", "longdash")) +
  theme(text=element_text(size=14)) + 
  geom_hline(yintercept=17/21, linetype="dashed", color = "grey") +
  theme(legend.title = element_blank()) +
  theme(legend.position = c(.70, .4))
ggsave(file.path("results", "unacknowledged_bias_grant.pdf"), height=4, width=4.5)



###Second, the augmented aldex simulation


fakeAldex.simulation <- function(d, n, seq.depth, pval = 0.05, prob = .9, test = "bayes_glm", n_samples = 2000){
  dd <- length(d)/2
  truth1 <- !(d[1:dd]==d[(dd+1):(2*dd)])##testing if the mean is different
  truth2 <- (1:dd)[truth1]##Locations of the differences
  dat <- create_true_abundances(d, n=n)
  rdat <- resample_data(dat, seq.depth=seq.depth)
  sample.totals=rowSums(dat[,-1])
  
  dd <- ncol(rdat)-1
  dfit <- run_deseq2(rdat)
  dfit <- summary_deseq2(dfit)
  dfit <- append_sig(dfit, function(x) sig_deseq2(x, pval=pval))
  
  afit <- run_aldex2(rdat)
  afit <- summary_aldex2(afit)
  afit <- append_sig(afit, function(x) sig_aldex2(x, pval=pval))
  
  upsilon = dd + 10
  Gamma = diag(c(5,1))
  Theta = matrix(0, dd-1, 2)
  Omega = diag(dd)
  G = cbind(diag(dd-1), -1)
  Xi = G%*%Omega%*%t(G)
  
  tfit.delta <- run_fakeAldex(rdat,
                         upsilon = upsilon, Gamma = Gamma,
                         Theta = Theta, Omega = Omega, Xi = Xi, n_samples = 2000, total_model = "logNormal", mean_lnorm = 0, sd_lnorm = 1e-3, test=test)
  tfit.naive <- run_fakeAldex(rdat,
                         upsilon = upsilon, Gamma = Gamma,
                         Theta = Theta, Omega = Omega, Xi = Xi, n_samples = 2000,  total_model = "logNormal", mean_lnorm = 0, sd_lnorm = 1, test=test)
  tfit.wk <- run_fakeAldex(rdat,
                      upsilon = upsilon, Gamma = Gamma,
                      Theta = Theta, Omega = Omega, Xi = Xi, n_samples = 2000,  total_model = "logNormal", mean_lnorm = miniclo(sample.totals), sd_lnorm = .5, sample.totals = sample.totals, test=test)
  tfit.coda <- run_fakeAldex(rdat,
                           upsilon = upsilon, Gamma = Gamma,
                           Theta = Theta, Omega = Omega, Xi = Xi, n_samples = 2000,  total_model = "logNormal", mean_lnorm = 0, sd_lnorm = 10, test=test)
  
  ##tfit.gm <- run_fakeAldex(rdat,
  ##                    upsilon = upsilon, Gamma = Gamma,
  ##                    Theta = Theta, Omega = Omega, Xi = Xi, n_samples = 2000, total_model = "gm", test=test)

  if(test == "bayes_glm"){
    tfit.delta <- summary_tram(tfit.delta, prob=prob)
    tfit.naive <- summary_tram(tfit.naive, prob=prob)
    tfit.wk <- summary_tram(tfit.wk, prob=prob)
    tfit.coda <- summary_tram(tfit.coda, prob=prob)
    tfit.delta <- append_sig(tfit.delta, sig_tram)
    tfit.naive <- append_sig(tfit.naive, sig_tram)
    tfit.wk <- append_sig(tfit.wk, sig_tram)
    tfit.coda <- append_sig(tfit.coda, sig_tram)
  } else{
    tfit.delta <- summary_aldex2(tfit.delta)
    tfit.naive <- summary_aldex2(tfit.naive)
    tfit.wk <- summary_aldex2(tfit.wk)
    tfit.coda <- summary_aldex2(tfit.coda)
    tfit.delta <- append_sig(tfit.delta, function(x) sig_aldex2(x, pval=pval))
    tfit.naive <- append_sig(tfit.naive, function(x) sig_aldex2(x, pval=pval))
    tfit.wk <- append_sig(tfit.wk, function(x) sig_aldex2(x, pval=pval))
    tfit.coda <- append_sig(tfit.coda, function(x) sig_aldex2(x, pval=pval))
  }
  
  rrs <- list(dat=rdat, dfit=dfit, afit=afit, tfit.delta = tfit.delta, tfit.naive = tfit.naive, tfit.wk = tfit.wk, tfit.coda = tfit.coda)
  ##rrs <- list(dat=rdat, dfit=dfit, afit=afit, tfit.naive = tfit.naive)
  
  p1 <- plot_count(dat)
  p2 <- plot_sig2(rrs, truth=truth2)
  p1 <- p1+theme(axis.title.x = element_blank(), 
                 axis.text.x = element_blank(), 
                 axis.ticks.x=element_blank(),
                 text = element_text(size=18))
  p <- plot_grid(p1, p2, nrow=2, align="v", rel_heights=c(1.7, 1))
  p
  #ggsave(p1, file.path("Results","aldex_deseq_failures.pdf"))
}# end of function



model.names <- c("dfit"="DESeq2",
                 "afit"="ALDEx2",
                 "tfit.delta" = "CLR",
                 "tfit.naive"= "Relaxed",
                 "tfit.wk" = "Informed",
                 "tfit.coda" = "CoDA")
model.name.levels <- c("CoDA","Informed", "Relaxed",  "CLR", "ALDEx2", "DESeq2")
# model.names <- c("dfit"="DESeq2",
#                  "afit"="ALDEx2",
#                  "tfit.naive"= "Relaxed")
# model.name.levels <- c("Relaxed", "ALDEx2", "DESeq2")


fakeAldex.simulation(d, n = 50, seq.depth = 5000, test = "bayes_glm", n_samples = 2000)


ggsave(file.path("results", "sim_matrixGraph.pdf"), height=4, width=7)
ggsave(file.path("results", "sim_matrixGraph_grant.pdf"), height=4, width=12)


plot_alpha <- function(d, n=50, seq.depth = 5000, alpha=seq(.01, 25, by=.5),
                       thresh=.9,...){
  dd <- length(d)/2
  truth1 <- !(d[1:dd]==d[(dd+1):(2*dd)])##testing if the mean is different
  truth2 <- (1:dd)[truth1]##Locations of the differences
  dat <- create_true_abundances(d, n=n)
  rdat <- resample_data(dat, seq.depth=seq.depth)
  sample.totals=rowSums(dat[,-1])
  
  dd <- ncol(rdat)-1
  
  coldata <- rdat[,"Condition",drop=F]
  countdata <- t(rdat[,-1,drop=F])
  rownames(coldata) <- colnames(countdata) <- paste0("n", 1:ncol(countdata))
  Y = as.matrix(countdata)
  X = as.character(coldata$Condition)
  
  alphaseq <- alpha
  fit <-run_fakeAldex(rdat,test = "t", total_model = "logNormal", mean_lnorm = 0, sd_lnorm = alpha[1])
  B <- matrix(NA, nrow = length(alpha), ncol = dd)
  pvals <- matrix(NA, nrow = length(alpha), ncol = dd)
  B[1,] <- fit$effect
  pvals[1,] <- fit$wi.ep
  if (length(alpha) > 1){
    for (i in 2:length(alpha)) {
      tmp = run_fakeAldex(rdat,test = "t", total_model = "logNormal", mean_lnorm = 0, sd_lnorm = alpha[i])
      B[i,] <- tmp$effect
      pvals[i,] <- tmp$wi.ep
    }
  }
  

  
  P = pvals %>% as.data.frame %>%
    as.data.frame() %>%
    mutate("alpha" = alpha) %>%
    dplyr::select(alpha, everything()) %>%
    pivot_longer(cols = !alpha, names_to = "Sequence", values_to = "pval")
  
  B %>% 
    as.data.frame() %>%
    mutate("alpha" = alpha) %>%
    dplyr::select(alpha, everything()) %>%
    pivot_longer(cols = !alpha, names_to = "Sequence", values_to = "Effect") %>%
    plyr::join(P, by = c("alpha", "Sequence")) %>%
    mutate("Sequence" = sub("V", "", Sequence)) %>%
    mutate("labl" = sub("V", "", Sequence)) %>%
    mutate("labl" = ifelse(labl %in% c(3, 4, 15, 21), labl, NA)) %>%
    ggplot(aes(x=alpha, y = Effect, group=Sequence)) +
    geom_line() +
    gghighlight((pval <= 0.05), use_direct_label  = FALSE) +
    gghighlight(!is.na(labl), unhighlighted_params = list(colour = NULL)) +
    geom_hline(yintercept=0, color="red", linetype = "dashed") +
    theme_bw() +
    ylab("Effect Size") +
    coord_cartesian(ylim = c(-3,2)) +
    scale_y_reverse() +
    xlab(TeX("$\\alpha$")) +
    theme(text = element_text(size=18))+
    theme(legend.position = "none") 
}
plot_alpha(d, alpha = seq(1e-3,5,by=.2))

ggsave(file.path("results", "sim_alphaGraph.pdf"), height=4, width=5)
ggsave(file.path("results", "sim_alphaGraph_grant.pdf"), height=4, width=7)

