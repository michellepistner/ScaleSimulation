###Code for the real data example

###Loading libraries
library(phyloseq)
library(fido)
library(tidyverse)
library(DESeq2)
library(ALDEx2)
library(gghighlight)
library(cowplot)
library(directlabels)
library(matrixNormal)
library(ggpattern)
library(LaplacesDemon)
library(ggsci)
library(mniw)


set.seed(1234)


# load data ---------------------------------------------------------------

# Load sequence count data
ps <- readRDS(file.path("data","phyloseq.rds"))

# Extract key elements as its hard to amalgamate (see next section) from
# w/in phyloseq
tax <- tax_table(ps)
otu <- as(t(otu_table(ps)), "matrix")

# Load flow cytometry data
flow <- read.csv(file=file.path("data", "counts_with_batch.csv"), header=TRUE) %>% 
  mutate(vessel = factor(vessel), 
         batch = factor(batch), 
         count = count *100000) # From jeff to get cell/ml



# preprocess data ---------------------------------------------------------

# Just take data from days 1 and 14 (only Hr4 from day 14) -- notice there are two sets from day 14
samples_to_keep <- str_subset(sample_names(ps), "^Day(1|14).Hr4\\.")
otu <- otu[,samples_to_keep]
flow_filtered <- filter(flow, day %in% c(1, 14))

# Filtering out taxa with < 2500 counts across all samples
# Putting these filtered taxa into an "other" category

otu %>% rowSums() %>% ecdf() %>% base::plot() %>% abline(v=2500)
to_amalgamate <- rowSums(otu) < 2500
tmp <- colSums(otu[to_amalgamate,])
otu <- rbind(otu[!to_amalgamate,], "other"=tmp)
otu_filtered <- otu

row.names(otu_filtered) = c(paste0("seq_", 1:(nrow(otu_filtered)-1)), "Other")

# Make a lookup table for the OTU data 
otu_filtered_metadata <- data.frame(sample = colnames(otu_filtered)) %>% 
  mutate(vessel = str_extract(sample, "(?<=V)[1-8]"),
         day = str_extract(sample, "(?<=^Day)[:digit:]+")) %>% 
  mutate(sample_id = paste0(day, "_", vessel)) 

# Model flow data ---------------------------------------------------------

flow_filtered <- flow_filtered %>% 
  mutate(sample_id = factor(paste0(flow_filtered$day,"_",flow_filtered$vessel)))


# Simple EDA: Relative abundance * flow -----------------------------------

flow %>%
  filter(day %in% c(1,14)) %>%
  mutate(day = factor(day)) %>% 
  ggplot(aes(x=day, y=count)) +
  geom_boxplot(fill="darkgrey") +
  theme_bw() +
  ylab("Bacterial Cells per mL") +
  xlab("Experiment Day")

ggsave(file.path("results", "vessel_1_flow_variation.pdf"), height=4, width=5)

# Subsetting the flow data 

flow_filtered_agg = flow_filtered %>%
  filter(vessel %in% c(3,4,5,6,7,8)) %>%
  dplyr::select(day, vessel, count, sample_id) %>%
  group_by(day, vessel, sample_id) %>%
  mutate(samp.var = sd(count)) %>%
  mutate(count = mean(count)) %>%
  ungroup() %>%
  unique()

# Set Priors  / data
Y <- otu_filtered
X <- model.matrix(~ base::I(vessel) + base::I(day) - 1,data = otu_filtered_metadata) %>% t()


options(ggrepel.max.overlaps = Inf)
plot_alpha <- function(Y, X, alpha=seq(0.01, 10, by=.5),
                       thresh=.9, upsilon = NULL, Theta = NULL, 
                       Gamma = NULL, Omega = NULL, Xi = NULL, total_model = "unif", sample.totals = NULL, sample = NULL, prob=.05, mean_lnorm = NULL, sd_lnorm = NULL, ...){
  dd <- nrow(Y)
  alphaseq <- alpha
  
  if(ncol(Y) != ncol(Y)){
    X = as.matrix(t(as.numeric(X)))}
  if(is.null(upsilon)){
    upsilon = dd + 3
  }
  if(is.null(Gamma)){
    Gamma = driver::bdiag(10*diag(6), 1)
  }
  if(is.null(Theta)){
    Theta = matrix(0, dd-1, nrow(X))
  }
  if(is.null(Omega)){
    Omega = diag(dd)
  }
  if(is.null(Xi)){
    G = cbind(diag(dd-1), -1)
    Xi = G%*%Omega%*%t(G)
  }
  
  ##Fitting and extracting results for the first alpha
  fit <-ssrv.mln(as.matrix(Y),X,covariate = "base::I(day)14",
                            upsilon = upsilon, Gamma = Gamma,
                            Theta = Theta, Omega = Omega, Xi = Xi, n_samples = 15000, sd_lnorm  = sqrt(alpha[1]^2/2),mean_lnorm = mean_lnorm, total_model = "logNormal", sample.totals =  sample.totals, sample = sample, prob = prob)

  B <- matrix(NA, nrow = length(alpha), ncol = dd)
  pvals <- matrix(NA, nrow = length(alpha), ncol = dd)
  B[1,] <- fit$mean/(fit$sd)
  pvals[1,] <- ifelse(sign(fit$low) == sign(fit$high), 1, 0)
  
  ##Repeating
  if (length(alpha) > 1){
    for (i in 2:length(alpha)) {
      tmp = ssrv.mln(as.matrix(Y),X,covariate = "base::I(day)14",
                                  upsilon = upsilon, Gamma = Gamma,
                                  Theta = Theta, Omega = Omega, Xi = Xi, n_samples = 15000, sd_lnorm  = sqrt(alpha[i]^2/2),mean_lnorm = mean_lnorm,  total_model = "logNormal", sample.totals =  sample.totals, sample = sample, prob = prob)
      B[i,] <- tmp$mean/(tmp$sd)
      pvals[i,] <- ifelse(sign(tmp$low) == sign(tmp$high), 1, 0)
      }
  }
  
  P = pvals %>% as.data.frame %>%
    as.data.frame() %>%
    mutate("alpha" = alpha) %>%
    dplyr::select(alpha, everything()) %>%
    pivot_longer(cols = !alpha, names_to = "Sequence", values_to = "pval")
  
  topValues <- colSums(pvals)[order(colSums(pvals), decreasing = TRUE)[1:5]]
  taxaToHighlight <- which(colSums(pvals) %in% topValues)
    
  p1 = B %>% 
    as.data.frame() %>%
    mutate("alpha" = alpha) %>%
    dplyr::select(alpha, everything()) %>%
    pivot_longer(cols = !alpha, names_to = "Sequence", values_to = "Effect") %>%
    plyr::join(P, by = c("alpha", "Sequence")) %>%
    mutate("Sequence" = sub("V", "", Sequence)) %>%
    mutate("labl" = sub("V", "", Sequence)) %>%
    mutate("labl" = ifelse(labl %in% taxaToHighlight, labl, NA)) %>%
    mutate(Effect = Effect) %>%
    ggplot(aes(x=alpha, y = Effect, group=Sequence)) +
    geom_line() +
    gghighlight((pval == TRUE), use_direct_label  = FALSE) +
    gghighlight(!is.na(labl), unhighlighted_params = list(colour = NULL)) +
    geom_hline(yintercept=0, color="red", linetype = "dashed") +
    geom_vline(xintercept = 0.5, color = "black", linetype = "dotted") +
    theme_bw() +
    ylab("Standardized LFC") +
    coord_cartesian(ylim = c(-4,2)) +
    scale_y_reverse() +
    xlab(expression(tau)) +
    theme(text = element_text(size=18))+
    theme(legend.position = "none") 
    
    p1
    return(list(p1 = p1))
}


plots = plot_alpha(Y,X, alpha=c(.1, .25, .5,1,  seq(2, 10, by=0.25)), total_model = "logNormal", mean_lnorm = rep(c(log(1),log(4)), each = 6), sample = 1:12, prob = .05)
plots$p1
ggsave(file.path("results", "realData_alpha.pdf"), height=4, width=7)


# Scale Simulation inference ----------------------------------------------------

# Set Priors  / data
Y <- otu_filtered
X <- model.matrix(~ base::I(vessel) + base::I(day)-1,data = otu_filtered_metadata) %>% t()
upsilon <- nrow(Y) + 3
Theta.t <- matrix(0, nrow(Y), nrow(X))
Theta.t[,7] <- 0
Gamma <- driver::bdiag(10*diag(6), 1)
Omega = diag(nrow(Y))
G = cbind(diag(nrow(Y)-1), -1)
Xi = G%*%Omega%*%t(G)

Theta <- G %*% Theta.t

fit_pim <- ssrv.mln(as.matrix(Y),X, covariate = rownames(X)[7],
                               upsilon = upsilon, Gamma = Gamma,
                               Theta = Theta, Omega = Omega, Xi = Xi, Theta.t = Theta.t, n_samples = 15000, total_model = "pim", prob = 0.05)


fit_design <- ssrv.mln(as.matrix(Y),X, covariate = rownames(X)[7],
                               upsilon = upsilon, Gamma = Gamma,
                               Theta = Theta, Omega = Omega, Xi = Xi, n_samples = 15000,  total_model = "logNormal",sd_lnorm = sqrt(.5^2/2), mean_lnorm=log(rep(c(1,4),each = 6)), sample = 1:12, prob = 0.05)


fit_gs <- ssrv.mln(as.matrix(Y),X, covariate = rownames(X)[7],
                               upsilon = upsilon, Gamma = Gamma,
                               Theta = Theta, Omega = Omega, Xi = Xi, n_samples = 15000, total_model = "flow", sample.totals =flow_filtered, sample = flow_filtered_agg$sample_id, prob = 0.05)


# Analysis using DESeq2 and Aldex2 ----------------------------------------

###Starting with DeSeq2

X.modelMat = t(X) %>% data.frame() %>% lapply(. %>% as.factor) %>% data.frame()

deseq.results = DESeqDataSetFromMatrix(countData=Y,
                                       colData=X.modelMat,
                                       design = ~ base..I.vessel.3 +base..I.vessel.4 +base..I.vessel.5+ base..I.vessel.6 + base..I.vessel.7 +base..I.day.14) ##Couldn't figure out how to automate this formula in a a way DeSeq would be happy
dds = DESeq(deseq.results, test="LRT", reduced =~   base..I.vessel.3 +base..I.vessel.4 +base..I.vessel.5+ base..I.vessel.6 + base..I.vessel.7 )
deseq_fit = results(dds)

deseq_results = deseq_fit %>%
  as.data.frame() %>% 
  rownames_to_column("category") %>% 
  dplyr::select(category, log2FoldChange, pvalue, lfcSE) %>% 
  mutate(low = log2FoldChange -1.96*lfcSE, 
         high = log2FoldChange + 1.96*lfcSE) %>% 
  mutate(mean=log2FoldChange) %>%
  filter(pvalue < .05)
deseq_results

###Aldex2 now
mmX <- model.matrix(~ base::I(vessel) + base::I(day),data = otu_filtered_metadata)
aldex_fit <- aldex(Y,mmX,"glm",mc.samples = 1000, denom = "all")

aldex_results = aldex_fit %>% 
  rownames_to_column("category") %>%
  dplyr::select(category,base..I.day.14.t.val, base..I.day.14.pval, base..I.day.14.pval.padj) %>%
  mutate(pval = base..I.day.14.pval) %>%
  mutate(padj = base..I.day.14.pval.padj) %>%
  filter(pval < .05)
aldex_results


# Summarise results and plot ----------------------------------------------

sig_gs <- fit_gs %>%
  mutate(sig = ifelse(sign(low) == sign(high), TRUE, FALSE)) %>%
  mutate(direction = -1*sign(mean)) %>% ## -1 needed to correct for direction changes to make comparable to aldex2/deseq2
  mutate(sig_code = ifelse(sig == TRUE, ifelse(direction ==1, "SigINC", "SigDEC"), "NS")) %>%
  mutate(res = ifelse(sig_code == "NS", "TN", "TP")) %>%
  dplyr::select(category, sig, sig_code,res) %>%
  mutate(Model = "SSRV - GS \n (Flow)     ")


sig_tram_design <- fit_design %>%
  mutate(sig = ifelse(sign(low) == sign(high), TRUE, FALSE)) %>%
  mutate(direction = -1*sign(mean)) %>% ## -1 needed to correct for direction changes to make comparable to aldex2/deseq2
  mutate(sig_code = ifelse(sig == TRUE, ifelse(direction ==1, "SigINC", "SigDEC"), "NS")) %>%
  mutate(res = "TN") %>%
  mutate(res = ifelse(((sig_code == sig_gs$sig_code) &(sig_gs$sig_code == "SigINC")) | ((sig_code == sig_gs$sig_code) &(sig_gs$sig_code == "SigDEC")), "TP", res)) %>%
  mutate(res = ifelse(((sig_code != "NS") &(sig_gs$sig_code == "NS")) | ((sig_code == "SigINC") &(sig_gs$sig_code == "SigDEC")) | ((sig_code == "SigDEC") &(sig_gs$sig_code == "SigINC")), "FP", res)) %>%
  mutate(res = ifelse(((sig_code == "NS") &(sig_gs$sig_code != "NS")), "FN", res)) %>%
  dplyr::select(category, sig, sig_code,res) %>%
  mutate(Model = "SSRV (Design) ")

sig_tram_pim <- fit_pim %>%
  mutate(sig = ifelse(sign(low) == sign(high), TRUE, FALSE)) %>%
  mutate(direction = -1*sign(mean)) %>% ## -1 needed to correct for direction changes to make comparable to aldex2/deseq2
  mutate(sig_code = ifelse(sig == TRUE, ifelse(direction ==1, "SigINC", "SigDEC"), "NS")) %>%
  mutate(res = "TN") %>%
  mutate(res = ifelse(((sig_code == sig_gs$sig_code) &(sig_gs$sig_code == "SigINC")) | ((sig_code == sig_gs$sig_code) &(sig_gs$sig_code == "SigDEC")), "TP", res)) %>%
  mutate(res = ifelse(((sig_code != "NS") &(sig_gs$sig_code == "NS")) | ((sig_code == "SigINC") &(sig_gs$sig_code == "SigDEC")) | ((sig_code == "SigDEC") &(sig_gs$sig_code == "SigINC")), "FP", res)) %>%
  mutate(res = ifelse(((sig_code == "NS") &(sig_gs$sig_code != "NS")), "FN", res)) %>%
  dplyr::select(category, sig, sig_code,res) %>%
  mutate(Model = "SSRV (PIM)")


deseq_results = deseq_fit %>%
  as.data.frame() %>%
  rownames_to_column("category") %>%
  mutate(sig = ifelse(padj < 0.05, TRUE, FALSE)) %>%
  mutate(direction = sign(log2FoldChange)) %>%
  mutate(sig_code = ifelse(sig == TRUE, ifelse(direction ==1, "SigINC", "SigDEC"), "NS")) %>%
  mutate(res = "TN") %>%
  mutate(res = ifelse(((sig_code == sig_gs$sig_code) &(sig_gs$sig_code == "SigINC")) | ((sig_code == sig_gs$sig_code) &(sig_gs$sig_code == "SigDEC")), "TP", res)) %>%
  mutate(res = ifelse(((sig_code != "NS") &(sig_gs$sig_code == "NS")) | ((sig_code == "SigINC") &(sig_gs$sig_code == "SigDEC")) | ((sig_code == "SigDEC") &(sig_gs$sig_code == "SigINC")), "FP", res)) %>%
  mutate(res = ifelse(((sig_code == "NS") &(sig_gs$sig_code != "NS")), "FN", res)) %>%
  dplyr::select(category, sig, sig_code,res) %>%
  mutate(Model = "DESeq2")


aldex_results = aldex_fit %>% 
  rownames_to_column("category") %>%
  mutate(sig = ifelse(base..I.day.14.pval.padj < 0.05, TRUE, FALSE)) %>%
  mutate(direction = sign(base..I.day.14.Est)) %>%
  mutate(sig_code = ifelse(sig == TRUE, ifelse(direction ==1, "SigINC", "SigDEC"), "NS")) %>%
  mutate(res = "TN") %>%
  mutate(res = ifelse(((sig_code == sig_gs$sig_code) &(sig_gs$sig_code == "SigINC")) | ((sig_code == sig_gs$sig_code) &(sig_gs$sig_code == "SigDEC")), "TP", res)) %>%
  mutate(res = ifelse(((sig_code != "NS") &(sig_gs$sig_code == "NS")) | ((sig_code == "SigINC") &(sig_gs$sig_code == "SigDEC")) | ((sig_code == "SigDEC") &(sig_gs$sig_code == "SigINC")), "FP", res)) %>%
  mutate(res = ifelse(((sig_code == "NS") &(sig_gs$sig_code != "NS")), "FN", res)) %>%
  dplyr::select(category, sig, sig_code,res) %>%
  mutate(Model = "ALDEx2")



##Generating the grid plot

sig.df <- rbind(aldex_results, deseq_results, sig_tram_design, sig_tram_pim, sig_gs) %>%
  filter(category != "Other")

sig.df$Model = factor(sig.df$Model, levels=c("ALDEx2", "DESeq2", "SSRV (Design) ", "SSRV (PIM)", "SSRV - GS \n (Flow)     "))
sig.df$Sequence = matrix(unlist(str_split(sig.df$category, "_")), ncol = 2, byrow = TRUE)[,2]

##No color labels
p2 = ggplot(sig.df, aes(x=Sequence, y=Model)) +
  geom_tile_pattern(aes(fill=res, pattern = res), color="darkgrey",pattern_fill = 'grey', pattern_colour  = 'grey', pattern_density = 0.015) +
  theme_minimal(18) +
  labs(title = "") +
  theme(panel.grid = element_blank(), 
        legend.title=element_blank(),
        legend.key.size = unit(.825, "cm")) +
  scale_pattern_manual(values = c(TP = "none", TN = "none", FP = "none", FN = "stripe")) +
  scale_fill_manual(values= c("white", "#fdae61", "white", "#2b83ba")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


ggsave(file.path("results", "model_comparison_flowData.pdf"), plot = p2, height=3, width=9)




####Repeated sampling
# Set Priors  / data
##num vessels
set.seed(2024)
num.vessels = c(6, 10,15,25,30,40,50)
replicates = 2

otu_filtered.day1 = otu_filtered[,1:6]
otu_filtered.day14 = otu_filtered[,7:12]

seq.names = rownames(Y)

## num false positives
fp = matrix(NA,nrow = length(num.vessels) * replicates, ncol = 6)
fp[,1] = rep(num.vessels, each = replicates)
colnames(fp) = c("num.vessels","ind", "Design","ALDEx2","DESeq2", "PIM")

## num true positives
tp = matrix(NA,nrow = length(num.vessels) * replicates, ncol = 6)
tp[,1] = rep(num.vessels, each = replicates)
colnames(tp) = c("num.vessels","ind", "Design","ALDEx2","DESeq2", "PIM")

## num true negatives
tn = matrix(NA,nrow = length(num.vessels) * replicates, ncol = 6)
tn[,1] = rep(num.vessels, each = replicates)
colnames(tn) = c("num.vessels","ind", "Design","ALDEx2","DESeq2", "PIM")

## num false negatives
fn = matrix(NA,nrow = length(num.vessels) * replicates, ncol = 6)
fn[,1] = rep(num.vessels, each = replicates)
colnames(fn) = c("num.vessels","ind", "Design","ALDEx2","DESeq2", "PIM")

ind <- 1
for(i in 1:length(num.vessels)){
  for(k in 1:replicates){
    # set the ind variable
    fn[ind,2] <- ind
    tn[ind,2] <- ind
    fp[ind,2] <- ind
    tp[ind,2] <- ind
      
    vessels.used = sample(3:8, num.vessels[i],replace = TRUE)

      
    # "-2" because the vessel numbers are 3-8
    otu.day1 = otu_filtered.day1[,c(vessels.used-2)]
    otu.day14 = otu_filtered.day14[,c(vessels.used-2)]
    
    flow_filtered_rep = c()
    for(j in 1:num.vessels[i]){  
      tmp = flow_filtered %>% filter(vessel == vessels.used[j]) %>%
        mutate(vessel = j) %>%
        mutate(sample_id = paste0(vessel, "_day", day))
      flow_filtered_rep = rbind(flow_filtered_rep, tmp)
    }
    colnames(otu.day1) = paste0("vessel_", 1:num.vessels[i], "_day1")
    colnames(otu.day14) = paste0("vessel_", 1:num.vessels[i], "_day14")
    
    Y = cbind(otu.day1,otu.day14)
    rownames(Y) = seq.names
      
    X.day1 = rbind(diag(num.vessels[i]), rep(0,num.vessels[i]))
    names(X.day1) = paste0("vessel_", 1:num.vessels[i], "_day1")
    X.day14 = rbind(diag(num.vessels[i]), rep(1,num.vessels[i]))
    names(X.day14) = paste0("vessel_", 1:num.vessels[i], "_day14")
    X = cbind(X.day1, X.day14)
    rownames(X) = c(paste0("base::I(vessel)", 1:num.vessels[i]), "base::I(day)14")
    
    upsilon <- nrow(Y) + 3
    Gamma <- driver::bdiag(10*diag(num.vessels[i]), 1) # Learn this function as well -- makes block diagonal matricies
    # This is the prior. Xi would technically be omega^comp
    Omega = diag(nrow(Y))
    Xi = G%*%Omega%*%t(G)
    Theta.t <- matrix(0, nrow(Y), nrow(X))
    Theta.t[,7] <- 0
    G = cbind(diag(nrow(Y)-1), -1)
    # prior for Theta^comp
    Theta <- G %*% Theta.t

    fit_gs <- ssrv.mln(as.matrix(Y),X, covariate = rownames(X)[nrow(X)],
                                  upsilon = upsilon, Gamma = Gamma,
                                  Theta = Theta, Omega = Omega, Xi = Xi, prob = 0.05, n_samples = 2000, total_model = "flow", sample.totals =flow_filtered_rep, sample = c(paste0( 1:num.vessels[i], "_day1"),paste0(1:num.vessels[i], "_day14")))
    sig_gs <- fit_gs %>%
      mutate(sig = ifelse(sign(low) == sign(high), TRUE, FALSE)) %>%
      mutate(direction = -1*sign(mean)) %>% ## -1 needed to correct for direction changes to make comparable to aldex2/deseq2
      mutate(sig_code = ifelse(sig == TRUE, ifelse(direction ==1, "SigINC", "SigDEC"), "NS"))
    
    fit_unif <- ssrv.mln(as.matrix(Y),X, covariate = rownames(X)[nrow(X)],
                                    upsilon = upsilon, Gamma = Gamma,
                                    Theta = Theta, Omega = Omega, Xi = Xi, prob = 0.05, n_samples = 2000,  total_model = "logNormal",sd_lnorm = sqrt(0.5^2/2), mean_lnorm= rep(c(log(1),log(4)), each = num.vessels[i]),  sample = 1:(2*num.vessels[i]))

    sig_unif <- fit_unif %>%
      mutate(sig = ifelse(sign(low) == sign(high), TRUE, FALSE)) %>%
      mutate(direction = -1*sign(mean)) %>% ## -1 needed to correct for direction changes to make comparable to aldex2/deseq2
      mutate(sig_code = ifelse(sig == TRUE, ifelse(direction ==1, "SigINC", "SigDEC"), "NS"))

    ##calculating statistics for design model
    tp[ind,3] <- sum((sig_unif$sig_code == sig_gs$sig_code) &(sig_gs$sig_code == "SigINC")) + sum((sig_unif$sig_code == sig_gs$sig_code) &(sig_gs$sig_code == "SigDEC"))
    fp[ind,3] <-  sum((sig_unif$sig_code != "NS") &(sig_gs$sig_code == "NS")) + sum((sig_unif$sig_code == "SigINC") &(sig_gs$sig_code == "SigDEC")) +  sum((sig_unif$sig_code == "SigDEC") &(sig_gs$sig_code == "SigINC"))
    fn[ind,3] <- sum((sig_unif$sig_code == "NS") &(sig_gs$sig_code != "NS"))
    tn[ind,3] <- sum((sig_unif$sig_code == "NS") &(sig_gs$sig_code == "NS"))
    
    # PIM
    fit_pim <- ssrv.mln(as.matrix(Y),X, covariate = rownames(X)[nrow(X)],
                        upsilon = upsilon, Gamma = Gamma,
                        Theta = Theta, Omega = Omega, Xi = Xi, Theta.t = Theta.t, n_samples = 2000, total_model = "pim", prob = 0.05)
    
    sig_pim <- fit_pim %>%
      mutate(sig = ifelse(sign(low) == sign(high), TRUE, FALSE)) %>%
      mutate(direction = -1*sign(mean)) %>% ## -1 needed to correct for direction changes to make comparable to aldex2/deseq2
      mutate(sig_code = ifelse(sig == TRUE, ifelse(direction ==1, "SigINC", "SigDEC"), "NS"))
    
    tp[ind,6] <- sum((sig_pim$sig_code == sig_gs$sig_code) &(sig_gs$sig_code == "SigINC")) + sum((sig_pim$sig_code == sig_gs$sig_code) &(sig_gs$sig_code == "SigDEC"))
    fp[ind,6] <-  sum((sig_pim$sig_code != "NS") &(sig_gs$sig_code == "NS")) + sum((sig_pim$sig_code == "SigINC") &(sig_gs$sig_code == "SigDEC")) +  sum((sig_pim$sig_code == "SigDEC") &(sig_gs$sig_code == "SigINC"))
    fn[ind,6] <- sum((sig_pim$sig_code == "NS") &(sig_gs$sig_code != "NS"))
    tn[ind,6] <- sum((sig_pim$sig_code == "NS") &(sig_gs$sig_code == "NS"))
    

    #running deseq2
    X.modelMat = t(X) %>% data.frame() %>% lapply(. %>% as.factor) %>% data.frame()
    
    deseq.design = as.formula(paste("~",paste(names(X.modelMat[,-1]), collapse = "+"),sep = " "))
    deseq.results = DESeqDataSetFromMatrix(countData=Y,
                                           colData=X.modelMat,
                                           design = deseq.design)
    deseq.design.red = as.formula(paste("~",paste(names(X.modelMat[,-c(1, ncol(X.modelMat))]), collapse = "+"),sep = " "))
    dds = DESeq(deseq.results, test="LRT", reduced = deseq.design.red)
    deseq_fit = results(dds)
    
    deseq_results = deseq_fit %>%
      as.data.frame() %>%
      mutate(sig = ifelse(padj < 0.05, TRUE, FALSE)) %>%
      mutate(direction = sign(log2FoldChange)) %>%
      mutate(sig_code = ifelse(sig == TRUE, ifelse(direction ==1, "SigINC", "SigDEC"), "NS"))
    
    tp[ind,5] <- sum((deseq_results$sig_code == sig_gs$sig_code) &(sig_gs$sig_code == "SigINC")) + sum((deseq_results$sig_code == sig_gs$sig_code) &(sig_gs$sig_code == "SigDEC"))
    fp[ind,5] <-  sum((deseq_results$sig_code != "NS") &(sig_gs$sig_code == "NS")) + sum((deseq_results$sig_code == "SigINC") &(sig_gs$sig_code == "SigDEC")) +  sum((deseq_results$sig_code == "SigDEC") &(sig_gs$sig_code == "SigINC"))
    fn[ind,5] <- sum((deseq_results$sig_code == "NS") &(sig_gs$sig_code != "NS"))
    tn[ind,5] <- sum((deseq_results$sig_code == "NS") &(sig_gs$sig_code == "NS"))
    
    
    ###Aldex2 now
    
    mmX <- t(X)
    aldex_fit <- aldex(Y,mmX,"glm",mc.samples = 128, denom = "all")
    
    aldex_results = aldex_fit %>% 
      rownames_to_column("category") %>%
      mutate(sig = ifelse(base..I.day.14.pval.padj < 0.05, TRUE, FALSE)) %>%
      mutate(direction = sign(base..I.day.14.Est)) %>%
      mutate(sig_code = ifelse(sig == TRUE, ifelse(direction ==1, "SigINC", "SigDEC"), "NS"))
    
    tp[ind,4] <- sum((aldex_results$sig_code == sig_gs$sig_code) &(sig_gs$sig_code == "SigINC")) + sum((aldex_results$sig_code == sig_gs$sig_code) &(sig_gs$sig_code == "SigDEC"))
    fp[ind,4] <-  sum((aldex_results$sig_code != "NS") &(sig_gs$sig_code == "NS")) + sum((aldex_results$sig_code == "SigINC") &(sig_gs$sig_code == "SigDEC")) +  sum((aldex_results$sig_code == "SigDEC") &(sig_gs$sig_code == "SigINC"))
    fn[ind,4] <- sum((aldex_results$sig_code == "NS") &(sig_gs$sig_code != "NS"))
    tn[ind,4] <- sum((aldex_results$sig_code == "NS") &(sig_gs$sig_code == "NS"))
    

    print(ind)
    ind <- ind+1
  }
  print(i)
}

tn <- tn %>%
  as.data.frame() %>%
  pivot_longer(!c(num.vessels, ind), names_to = "method", values_to = "tn")

typei <- fp %>%
  as.data.frame() %>%
  pivot_longer(!c(num.vessels, ind), names_to = "method", values_to = "fp") %>%
  plyr::join(tn, by =c("num.vessels", "ind", "method")) %>%
  mutate(fdr = ifelse((fp+tn) > 0, fp/(fp+tn), 0)) %>%
  dplyr::select(-ind) %>%
  group_by(num.vessels, method) %>%
  mutate(mean = mean(fdr, na.rm=TRUE)) %>%
  mutate(meanfp = mean(fp, na.rm= TRUE)) %>%
  mutate(sd = sd(fdr, na.rm = TRUE)) %>%  
  mutate(sdfp = sd(fp, na.rm = TRUE)) %>%
  ungroup()  %>%
  dplyr::select(-c(fp, tn, fdr)) %>%
  unique()

##plotting
plot1 <- ggplot(typei, aes(x = num.vessels, y = mean, group = method,linetype = method, color = method, fill = method)) +
  #geom_ribbon(aes(ymin = mean-sd, ymax=mean+sd), alpha = .1, colour = NA) +
  scale_color_manual(values=c("#AF4BCE", "#EA7369", "#2b83ba", "darkgreen")) + 
  scale_linetype_manual(values = c("dotted", "twodash", "longdash", "solid")) +
  geom_line(alpha = 1, lwd = 1.1) +
  theme_bw() +
  xlab("Number of Vessels") +
  ylab("Type-I Error") +
  theme(legend.title = element_blank(),
        text = element_text(size=16))
plot1

ggsave(file.path("results", "typei-by-method-gut-data.pdf"), height=4, width=4.5)

plot1 <- ggplot(typei, aes(x = num.vessels, y = mean, group = method,linetype = method, color = method, fill = method)) +
  geom_ribbon(aes(ymin = mean-sd, ymax=mean+sd), alpha = .1, colour = NA) +
  scale_color_manual(values=c("#AF4BCE", "#EA7369", "#2b83ba", "darkgreen")) + 
  scale_fill_manual(values=c("#AF4BCE", "#EA7369", "#2b83ba", "darkgreen")) + 
  scale_linetype_manual(values = c("dotted", "twodash", "longdash", "solid")) +
  geom_line(alpha = 1, lwd = 1.1) +
  theme_bw() +
  xlab("Number of Vessels") +
  ylab("Type-I Error") +
  theme(legend.title = element_blank(),
        text = element_text(size=16))
plot1

ggsave(file.path("results", "typei-by-method-gut-data_errors.pdf"), height=4, width=4.5)


plot2 <- ggplot(typei, aes(x = num.vessels, y = meanfp, group = method, linetype = method, color = method, fill = method)) +
  #geom_ribbon(aes(ymin = meanfp-sdfp, ymax=meanfp+sdfp), alpha = .1, colour = NA) +
  scale_color_manual(values=c("#AF4BCE", "#EA7369", "#2b83ba", "darkgreen")) + 
  scale_linetype_manual(values = c("dotted", "twodash", "longdash", "solid")) +
  geom_line(alpha = 1, lwd = 1.1) +
  theme_bw() +
  xlab("Number of Vessels") +
  ylab("Number of False Positives") +
  theme(legend.title = element_blank(),
        text = element_text(size=16))
plot2

ggsave(file.path("results", "fp-by-method-gut-data.pdf"), height=4, width=4.5)

plot2 <- ggplot(typei, aes(x = num.vessels, y = meanfp, group = method, linetype = method, color = method, fill = method)) +
  geom_ribbon(aes(ymin = meanfp-sdfp, ymax=meanfp+sdfp), alpha = .1, colour = NA) +
  scale_color_manual(values=c("#AF4BCE", "#EA7369", "#2b83ba", "darkgreen")) + 
  scale_fill_manual(values=c("#AF4BCE", "#EA7369", "#2b83ba", "darkgreen")) + 
  scale_linetype_manual(values = c("dotted", "twodash", "longdash", "solid")) +
  geom_line(alpha = 1, lwd = 1.1) +
  theme_bw() +
  xlab("Number of Vessels") +
  ylab("Number of False Positives") +
  theme(legend.title = element_blank(),
        text = element_text(size=16))
plot2

ggsave(file.path("results", "fp-by-method-gut-data-errors.pdf"), height=4, width=4.5)


tp <- tp %>%
  as.data.frame() %>%
  pivot_longer(!c(num.vessels, ind), names_to = "method", values_to = "tp")

fdr <- fp %>%
  as.data.frame() %>%
  pivot_longer(!c(num.vessels, ind), names_to = "method", values_to = "fp") %>%
  plyr::join(tp, by =c("num.vessels", "ind", "method")) %>%
  mutate(fdr = ifelse((fp+tp) > 0, fp/(fp+tp), 0)) %>%
  dplyr::select(-ind) %>%
  group_by(num.vessels, method) %>%
  mutate(mean = mean(fdr, na.rm = TRUE)) %>%
  mutate(sd = sd(fdr, na.rm = TRUE)) %>%
  ungroup() %>%
  dplyr::select(-c(fp, tp, fdr)) %>%
  unique()

##plotting
plot3 <- ggplot(fdr, aes(x = num.vessels, y = mean, group = method, linetype = method, color = method, fill = method)) +
  #geom_ribbon(aes(ymin = mean-sd, ymax=mean+sd), alpha = .1, colour = NA) +
  scale_color_manual(values=c("#AF4BCE", "#EA7369", "#2b83ba", "darkgreen")) + 
  scale_linetype_manual(values = c("dotted", "twodash", "longdash", "solid")) +
  geom_line(alpha = 1, lwd = 1.1) +
  theme_bw() +
  xlab("Number of Vessels") +
  ylab("False Discovery Rate") +
  theme(legend.title = element_blank(),
        text = element_text(size=16))
plot3

ggsave(file.path("results", "fdr-by-method-gut-data.pdf"), height=4, width=4.5)

plot3 <- ggplot(fdr, aes(x = num.vessels, y = mean, group = method, linetype = method, color = method, fill = method)) +
  geom_ribbon(aes(ymin = mean-sd, ymax=mean+sd), alpha = .1, colour = NA) +
  scale_color_manual(values=c("#AF4BCE", "#EA7369", "#2b83ba", "darkgreen")) + 
  scale_fill_manual(values=c("#AF4BCE", "#EA7369", "#2b83ba", "darkgreen")) + 
  scale_linetype_manual(values = c("dotted", "twodash", "longdash", "solid")) +
  geom_line(alpha = 1, lwd = 1.1) +
  theme_bw() +
  xlab("Number of Vessels") +
  ylab("False Discovery Rate") +
  theme(legend.title = element_blank(),
        text = element_text(size=16))
plot3

ggsave(file.path("results", "fdr-by-method-gut-data-errors.pdf"), height=4, width=4.5)

# Type II error
typeii <- fn %>%
  as.data.frame() %>%
  pivot_longer(!c(num.vessels, ind), names_to = "method", values_to = "fn") %>%
  plyr::join(tp, by =c("num.vessels", "ind", "method")) %>%
  mutate(typeii = ifelse((fn+tp) > 0, fn/(fn+tp), 0)) %>%
  dplyr::select(-ind) %>%
  group_by(num.vessels, method) %>%
  mutate(mean = mean(typeii, na.rm = TRUE)) %>%
  mutate(sd = sd(typeii, na.rm = TRUE)) %>%
  ungroup() %>%
  dplyr::select(-c(fn, tp, typeii)) %>%
  unique()
  

##plotting##plottingTRUE
plot4 <- ggplot(typeii, aes(x = num.vessels, y = mean, linetype = method, group = method, color = method, fill = method)) +
  #geom_ribbon(aes(ymin = mean-sd, ymax=mean+sd), alpha = .1, colour = NA) +
  scale_color_manual(values=c("#AF4BCE", "#EA7369", "#2b83ba", "darkgreen")) + 
  scale_linetype_manual(values = c("dotted", "twodash", "longdash", "solid")) +
  geom_line(alpha = 1, lwd = 1.1) +
  theme_bw() +
  xlab("Number of Vessels") +
  ylab("Type-II Error") +
  theme(legend.title = element_blank(),
        text = element_text(size=16))
plot4

ggsave(file.path("results", "typeii-by-method-gut-data.pdf"), height=4, width=4.5)

plot4 <- ggplot(typeii, aes(x = num.vessels, y = mean, linetype = method, group = method, color = method, fill = method)) +
  geom_ribbon(aes(ymin = mean-sd, ymax=mean+sd), alpha = .1, colour = NA) +
  scale_color_manual(values=c("#AF4BCE", "#EA7369", "#2b83ba", "darkgreen")) + 
  scale_linetype_manual(values = c("dotted", "twodash", "longdash")) +
  scale_fill_manual(values=c("#AF4BCE", "#EA7369", "#2b83ba", "darkgreen")) + 
  geom_line(alpha = 1, lwd = 1.1) +
  theme_bw() +
  xlab("Number of Vessels") +
  ylab("Type-II Error") +
  theme(legend.title = element_blank(),
        text = element_text(size=16))
plot4

ggsave(file.path("results", "typeii-by-method-gut-data-errors.pdf"), height=4, width=4.5)
