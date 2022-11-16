##Code for the simulation example

###Loading libraries
library(tidyverse)
library(DESeq2)
library(ALDEx2)
library(LaplacesDemon)
library(ggpattern)
library(cowplot)
library(gghighlight)
library(directlabels)
library(latex2exp)

set.seed(1234)

###Setting the data parameters for all simulations
d <- c(4000, 4000, 4000, 4000, 4000, 400,400,400,400,4000,400,500,500,500,400,400,400,400,400,400, # Pre
       4000, 4000, 3000, 2000, 4000, 400,400,400,400,4000,400,500,500,500,200,400,400,400,400,100) # Post
###Half are for the case, half for control
dd = length(d)/2

##Finding which are truly DE
truth1 <- !(d[1:dd]==d[(dd+1):(2*dd)])##testing if the mean is different


##Helper function for the ALDEx2 + DESeq2 analysis of FDR
aldexDeseq_analysis <- function(d, n, seq.depth, pval = 0.05, prob = .99){
  ##Finding the number of taxa and which are different
  dd <- length(d)/2
  truth1 <- !(d[1:dd]==d[(dd+1):(2*dd)])##testing if the mean is different

  ##Creating the true data
  dat <- create_true_abundances(d, n=n)
  ##Resampling the data
  rdat <- resample_data(dat, seq.depth=seq.depth)
  
  ##DESeq2 analysis: run on the resampled data. Then, find those that are significant
  dfit <- run_deseq2(rdat)
  dfit <- summary_deseq2(dfit)
  dfit <- append_sig(dfit, function(x) sig_deseq2(x, pval=pval))
  
  ##ALDEx2 analysis: run on the resampled data. Then, find those that are signficant
  afit <- run_aldex2(rdat)
  afit <- summary_aldex2(afit)
  afit <- append_sig(afit, function(x) sig_aldex2(x, pval=pval))
  
  ##Find the true positives and false positives for each
  fp.deseq = sum(dfit$sig == TRUE & truth1 == FALSE)
  tp.deseq = sum(dfit$sig == TRUE & truth1 == TRUE)
  
  fp.aldex = sum(afit$sig == TRUE & truth1 == FALSE)
  tp.aldex = sum(afit$sig == TRUE & truth1 == TRUE)
  
  
  return(list("fp_deseq" = fp.deseq, "tp_deseq" = tp.deseq, "fp_aldex" = fp.aldex, "tp_aldex" = tp.aldex))
  
}# end of function

##Running the above helper function
model.names <- c("dfit"="DESeq2",
                 "afit"="ALDEx2")
model.name.levels <- c("DESeq2", "ALDEx2")

##Number of conditions to run it on
vals = c(5,10,25,50,75,100, 125, 150, 200)

fdr.aldex = rep(NA, length(vals))
fdr.deseq = rep(NA, length(vals))

##Repeat for each sample size
for(i in 1:length(vals)){
  ##Run the analysis
  res = aldexDeseq_analysis(d, vals[i], 5000)
  
  ##Calculate FDR
  fdr.aldex[i] = res$fp_aldex/(res$fp_aldex + res$tp_aldex)
  fdr.deseq[i] = res$fp_deseq/(res$fp_deseq + res$tp_deseq)

  print(i)
}

##Now, graphing the FDR rates by methods...
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

##Second, showing how ALDEx2 matches our SSRV
aldexSSRV_analysis <- function(d, n, seq.depth, pval = 0.05, prob = .99){
  ##Find the number of samples per condition
  dd <- length(d)/2
  truth1 <- !(d[1:dd]==d[(dd+1):(2*dd)])##testing if the mean is different

  ##Create the true abundance data set; then create the resampled data set
  dat <- create_true_abundances(d, n=n)
  rdat <- resample_data(dat, seq.depth=seq.depth)
  
  ##Running ALDEx2
  afit <- run_aldex2(rdat)
  afit <- summary_aldex2(afit)
  afit <- append_sig(afit, function(x) sig_aldex2(x, pval=pval))
  
  ##Repeating with the SSRV version
  ssrv <- run_fakeAldex(rdat, n_samples = 1000, gamma = 1e-3)
  ssrv <- summary_aldex2(ssrv)
  ssrv <- append_sig(ssrv, function(x) sig_aldex2(x, pval = pval))
  
  rrs <- list(dat=rdat, ssrv=ssrv, afit=afit)
  
  ##Finding the number of false positives and true positives for each
  fp.ssrv = sum(ssrv$sig == TRUE & truth1 == FALSE)
  tp.ssrv = sum(ssrv$sig == TRUE & truth1 == TRUE)
  
  fp.aldex = sum(afit$sig == TRUE & truth1 == FALSE)
  tp.aldex = sum(afit$sig == TRUE & truth1 == TRUE)
  
  return(list("fp_ssrv" = fp.ssrv, "tp_ssrv" = tp.ssrv, "pval_aldex" = afit$padj, "fp_aldex" = fp.aldex, "tp_aldex" = tp.aldex))
  
}# end of function

##Running now for the SSRV + ALDEx2
model.names <- c("ssrv"="SSRV (Geo. Mean)",
                 "afit"="ALDEx2")
model.name.levels <- c("SSRV (Geo. Mean)", "ALDEx2")


vals = c(5,10,25,50,75,100, 125, 150, 200)

fdr.aldex = rep(NA, length(vals))
fdr.ssrv = rep(NA, length(vals))

##Repeating for every sample size above
for(i in 1:length(vals)){
  res = aldexSSRV_analysis(d, vals[i], 5000)
  
  ##Calculating FDR rates
  fdr.aldex[i] = res$fp_aldex/(res$fp_aldex + res$tp_aldex)
  fdr.ssrv[i] = res$fp_ssrv/(res$fp_ssrv + res$tp_ssrv)
  
  print(i)
}

##Graphing and saving the graph
fdr.all = data.frame(vals = rep(vals,2), fdr = c(c(fdr.aldex),c(fdr.ssrv)), method = rep(c("ALDEx2", "SSRV (Geo. Mean)"), each = length(vals)))
fdr.all$method = as.factor(fdr.all$method)

ggplot(fdr.all, aes(x=vals, y=fdr, color=method, fill = method, linetype = method)) +
  geom_line(alpha = 1, lwd = 1.1, position = position_dodge(width = 0.5)) +
  theme_bw()+
  xlim(c(0,200)) +
  xlab("Sample Size") +
  ylab("False Discovery Rate") + 
  scale_color_manual(values=c("#e66101", "#5e3c99")) + 
  scale_linetype_manual(values = c("twodash", "dotted")) +
  theme(text=element_text(size=14)) + 
  geom_hline(yintercept=16/20, linetype="dashed", color = "grey") +
  theme(legend.title = element_blank())
ggsave(file.path("results", "recreating_aldex2.pdf"), height=3, width=6)


###Third, the SSRV simulation over different total models
fakeAldex.simulation <- function(d, n, seq.depth, pval = 0.05, prob = .9, test = "t", n_samples = 2000){
  dd <- length(d)/2
  truth1 <- !(d[1:dd]==d[(dd+1):(2*dd)])##testing if the mean is different
  dat <- create_true_abundances(d, n=n)
  rdat <- resample_data(dat, seq.depth=seq.depth)

  ##DESeq2
  dfit <- run_deseq2(rdat)
  dfit <- summary_deseq2(dfit)
  dfit <- append_sig(dfit, function(x) sig_deseq2(x, pval=pval))
  
  ##ALDEx2
  afit <- run_aldex2(rdat)
  afit <- summary_aldex2(afit)
  afit <- append_sig(afit, function(x) sig_aldex2(x, pval=pval))
  
  ##SSRV models
  ssrv.delta <- run_fakeAldex(rdat, n_samples = 2000, gamma = 1e-3)
  ssrv.delta <- summary_aldex2(ssrv.delta)
  ssrv.delta <- append_sig(ssrv.delta, function(x) sig_aldex2(x, pval = pval))
  
  ssrv.noise <- run_fakeAldex(rdat, n_samples = 2000, gamma = .5)
  ssrv.noise <- summary_aldex2(ssrv.noise)
  ssrv.noise <- append_sig(ssrv.noise, function(x) sig_aldex2(x, pval = pval))
  
  ssrv.coda <- run_fakeAldex(rdat, n_samples = 2000, gamma = 10)
  ssrv.coda <- summary_aldex2(ssrv.coda)
  ssrv.coda <- append_sig(ssrv.coda, function(x) sig_aldex2(x, pval = pval))
  
  rrs <- list(dat=rdat, dfit = dfit, afit = afit, ssrv.delta = ssrv.delta, ssrv.noise = ssrv.noise, ssrv.coda = ssrv.coda)
  p1 <- plot_count(dat)
  p2 <- plot_sig2(rrs, truth=truth2)
  p1 <- p1+theme(axis.title.x = element_blank(), 
                 axis.text.x = element_blank(), 
                 axis.ticks.x=element_blank(),
                 text = element_text(size=14))
  p <- plot_grid(p1, p2, nrow=2, align="v", rel_heights=c(1.7, 1))
  p
}# end of function


##Models graphed
model.names <- c("dfit"="DESeq2",
                 "afit"="ALDEx2",
                 "ssrv.delta" = "CLR",
                 "ssrv.noise"= "Relaxed",
                 "ssrv.coda" = "CoDA")
model.name.levels <- c("CoDA", "Relaxed",  "CLR", "ALDEx2", "DESeq2")

##Running the simulation. Output is a graph
fakeAldex.simulation(d, n = 50, seq.depth = 5000, test = "t", n_samples = 2000)

##Saving the graph
ggsave(file.path("results", "sim_matrixGraph.pdf"), height=4, width=4.5)


plot_alpha <- function(d, n=50, seq.depth = 5000, alpha=seq(.01, 25, by=.5),
                       thresh=.9,...){
  dd <- length(d)/2
  truth1 <- !(d[1:dd]==d[(dd+1):(2*dd)])##testing if the mean is different
  dat <- create_true_abundances(d, n=n)
  rdat <- resample_data(dat, seq.depth=seq.depth)

  coldata <- rdat[,"Condition",drop=F]
  countdata <- t(rdat[,-1,drop=F])
  rownames(coldata) <- colnames(countdata) <- paste0("n", 1:ncol(countdata))
  Y = as.matrix(countdata)
  X = as.character(coldata$Condition)
  
  ##Running for the first alpha and extracting the effect and p-value
  fit <-run_fakeAldex(rdat, n_samples = 2000, gamma = alpha[1])
  B <- matrix(NA, nrow = length(alpha), ncol = dd)
  pvals <- matrix(NA, nrow = length(alpha), ncol = dd)
  B[1,] <- fit$effect ##Creating and effect size: mean adjusted by length of CI
  pvals[1,] <- ifelse(fit$we.eBH < 0.05, TRUE, FALSE)
  
  ##Now repeating for the rest of the values of alpha
  if (length(alpha) > 1){
    for (i in 2:length(alpha)) {
      tmp = run_fakeAldex(rdat, n_samples = 2000, gamma = alpha[i])
      B[i,] <- tmp$effect
      pvals[i,] <- ifelse(tmp$we.eBH < 0.05, TRUE, FALSE)
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
    mutate("labl" = ifelse(labl %in% c(3, 4, 15, 20), labl, NA)) %>%
    ggplot(aes(x=alpha, y = Effect, group=Sequence)) +
    geom_line() +
    gghighlight((pval == TRUE), use_direct_label  = FALSE) +
    gghighlight(!is.na(labl), unhighlighted_params = list(colour = NULL)) +
    geom_hline(yintercept=0, color="red", linetype = "dashed") +
    theme_bw() +
    ylab("Effect Size") +
    coord_cartesian(ylim = c(-10,6)) +
    scale_y_reverse() +
    xlab(TeX("$\\alpha$")) +
    theme(text = element_text(size=18))+
    theme(legend.position = "none") 
}

##Running, plotting, and saving
plot_alpha(d, alpha = seq(1e-3,4,by=.1))

ggsave(file.path("results", "sim_alphaGraph.pdf"), height=4, width=5)
