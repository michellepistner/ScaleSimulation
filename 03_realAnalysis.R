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
library(latex2exp)

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

# filter samples
otu %>% colSums() %>% ecdf() %>% base::plot(xlim=c(0,50000))
otu <- otu[,colSums(otu) > 10000] 

# Just take data from days 1 and 14 (only Hr4 from day 14) -- notice there are two sets from day 14
samples_to_keep <- str_subset(sample_names(ps), "^Day(1|14).Hr4\\.")
otu <- otu[,samples_to_keep]
flow_filtered <- filter(flow, day %in% c(1, 14))

(rowSums(otu > 3)/ncol(otu)) %>% ecdf  %>% base::plot(xlim=c(0,0.4))
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

# Set Priors  / data
Y <- otu_filtered
X <- model.matrix(~ base::I(vessel) + base::I(day) - 1,data = otu_filtered_metadata) %>% t()


options(ggrepel.max.overlaps = Inf)
plot_alpha <- function(Y, X, alpha=seq(0.01, 10, by=.5),
                       thresh=.9, upsilon = NULL, Theta = NULL, 
                       Gamma = NULL, Omega = NULL, Xi = NULL, total_model = "unif", sample.totals = NULL, sample = NULL, prob=.975, mean_lnorm = NULL, sd_lnorm = NULL, ...){
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
                            Theta = Theta, Omega = Omega, Xi = Xi, n_samples = 15000, sd_lnorm  = sqrt(alpha[1]^2/2),mean_lnorm = mean_lnorm, total_model = "logNormal", sample.totals =  sample.totals, sample = sample)

  B <- matrix(NA, nrow = length(alpha), ncol = dd)
  pvals <- matrix(NA, nrow = length(alpha), ncol = dd)
  B[1,] <- fit$mean/(fit$high - fit$low)
  pvals[1,] <- ifelse(sign(fit$low) == sign(fit$high), 1, 0)
  
  ##Repeating
  if (length(alpha) > 1){
    for (i in 2:length(alpha)) {
      tmp = ssrv.mln(as.matrix(Y),X,covariate = "base::I(day)14",
                                  upsilon = upsilon, Gamma = Gamma,
                                  Theta = Theta, Omega = Omega, Xi = Xi, n_samples = 15000, sd_lnorm  = sqrt(alpha[i]^2/2),mean_lnorm = mean_lnorm,  total_model = "logNormal", sample.totals =  sample.totals, sample = sample)
      B[i,] <- tmp$mean/(tmp$high - tmp$low)
      pvals[i,] <- ifelse(sign(tmp$low) == sign(tmp$high), 1, 0)
      }
  }
  
  P = pvals %>% as.data.frame %>%
    as.data.frame() %>%
    mutate("alpha" = alpha) %>%
    dplyr::select(alpha, everything()) %>%
    pivot_longer(cols = !alpha, names_to = "Sequence", values_to = "pval")
  
  topValues <- colSums(pvals)[order(colSums(pvals), decreasing = TRUE)[1:10]]
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
    theme_bw() +
    ylab("Effect Size") +
    coord_cartesian(ylim = c(-10,6)) +
    scale_y_reverse() +
    xlab(TeX("$\\alpha$")) +
    theme(text = element_text(size=18))+
    theme(legend.position = "none") 
    
    p1
    return(list(p1 = p1))
}

flow_filtered_agg = flow_filtered %>%
  filter(vessel %in% c(3,4,5,6,7,8)) %>%
  dplyr::select(day, vessel, count, sample_id) %>%
  group_by(day, vessel, sample_id) %>%
  mutate(samp.var = sd(count)) %>%
  mutate(count = mean(count)) %>%
  ungroup() %>%
  unique()

plots = plot_alpha(Y,X, alpha=c(.1, .25, .5,1,  seq(2, 10, by=0.25)), total_model = "logNormal", mean_lnorm = rep(c(log(1),log(4)), each = 6), sample = 1:12, prob = .975)
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
  dplyr::select(category, log2FoldChange, padj, lfcSE) %>% 
  mutate(low = log2FoldChange -1.96*lfcSE, 
         high = log2FoldChange + 1.96*lfcSE) %>% 
  mutate(mean=log2FoldChange) %>%
  filter(padj < .05)
deseq_results

###Aldex2 now
mmX <- model.matrix(~ base::I(vessel) + base::I(day),data = otu_filtered_metadata)
aldex_fit <- aldex(Y,mmX,"glm",mc.samples = 1000, denom = "all")

aldex_results = aldex_fit %>% 
  rownames_to_column("category") %>%
  dplyr::select(category,model.base..I.day.14.t.value, model.base..I.day.14.Pr...t.., model.base..I.day.14.Pr...t...BH) %>%
  mutate(pval = model.base..I.day.14.Pr...t..) %>%
  mutate(padj = model.base..I.day.14.Pr...t...BH) %>%
  filter(pval < .05)
aldex_results


# Summarise results and plot ----------------------------------------------

sig_tram_design = sig_tram(fit_design)

sig_tram_gs = sig_tram(fit_gs)

sig_tram_pim = sig_tram(fit_pim)

###Extract anything that's significant anywhere across all models
sig.values = c(aldex_results$category,deseq_results$category, sig_tram_design$category, sig_tram_pim$category, sig_tram_gs$category) %>% unique

###Some light processing to make it more useful

truth.pos = sig_tram_gs$category
truth.pos = sub("seq_", "", truth.pos)

##Generating the grid plot
q=length(sig.values)

sig.df = data.frame("Sequence" = rep(str_split_fixed(sig.values, "\\_", 2)[,2],5))
sig.df = sig.df %>%
  mutate(Sequence = ifelse(Sequence == "", "Other", Sequence)) %>%
  mutate(true.pos = ifelse(Sequence %in% truth.pos, 1, 0)) %>%
  mutate(Model = c(rep("ALDEx2", q),rep("DESeq2", q), rep("SSRV (Design) ", q),  rep("SSRV (PIM)", q), rep("SSRV - GS \n (Flow)     ", q))) %>%
  mutate(sigcode = c(ifelse(sig.values %in% aldex_results$category, 1, 0),ifelse( sig.values %in% deseq_results$category, 1, 0),ifelse( sig.values %in% sig_tram_design$category, 1, 0), ifelse( sig.values %in% sig_tram_pim$category, 1, 0), ifelse( sig.values %in% sig_tram_gs$category, 1, 0))) %>%
  mutate(res = ifelse(true.pos == 1 & sigcode == 1, "TP", NA)) %>%
  mutate(res = ifelse(true.pos == 0 & sigcode == 0, "TN", res)) %>%
  mutate(res = ifelse(true.pos == 1 & sigcode == 0, "FN", res)) %>%
  mutate(res = ifelse(true.pos == 0 & sigcode == 1, "FP", res)) %>%
  mutate(sigcode = factor(sigcode, levels = list("Non. Sig."="0", "Sig."="1"))) %>%
  mutate(Sequence = as.numeric(sig.df$Sequence)) %>%
  arrange(Sequence) %>%
  mutate(Sequence = as.factor(Sequence)) %>%
  filter(!is.na(Sequence))


sig.df$Model = factor(sig.df$Model, levels=c("ALDEx2", "DESeq2", "SSRV (CLR) ", "SSRV (CoDA) ", "SSRV (Design) ", "SSRV (PIM)", "SSRV - GS \n (Flow)     "))
sig.df$Sequence = as.character(sig.df$Sequence)
sig.df$Sequence = factor(sig.df$Sequence, levels =c(as.character(sort(as.numeric(str_split_fixed(sig.values, "\\_", 2)[,2]))), "Other"))

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
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), text = element_text(size=16))

p2


ggsave(file.path("results", "model_comparison_flowData.pdf"), height=4, width=5)




####Repeated sampling
# Set Priors  / data
##num vessels

num.vessels = c(6, 10,15,25,30,40,50)

size.day1 = colSums(otu_filtered[,1:6])
size.day14 = colSums(otu_filtered[,7:12])
  
otu_filtered.day1 = otu_filtered[,1:6]
otu_filtered.day14 = otu_filtered[,7:12]

seq.names = rownames(Y)
fdr = matrix(NA,nrow = length(num.vessels), ncol = 4)
fdr[,1] = num.vessels
names(fdr) = c("num.vessels","relaxed","aldex2","deseq2")
replicates = 10
for(i in 1:length(num.vessels)){
  tmp.fdr = rep(0,3)
  for(k in 1:replicates){
    if(num.vessels[i] != 6){
      vessels.used = sample(1:6, num.vessels[i],replace = TRUE)
    } else{
      vessels.used = sample(1:6, num.vessels[i],replace = FALSE)
    }
    otu.day1 = matrix(nrow = length(seq.names))
    otu.day14 = matrix(nrow = length(seq.names))
    flow_filtered_rep = c()
    for(j in 1:num.vessels[i]){
      otu.day1 = cbind(otu.day1, c(rmultinom(1,size.day1[vessels.used[j]], otu_filtered.day1[,vessels.used[j]])))
      otu.day14 = cbind(otu.day14, c(rmultinom(1,size.day14[vessels.used[j]], otu_filtered.day14[,vessels.used[j]])))
      
      tmp = flow_filtered %>% filter(vessel == vessels.used[j]) %>%
        mutate(vessel = j) %>%
        mutate(sample_id = paste0(vessel, "_day", day))
      flow_filtered_rep = rbind(flow_filtered_rep, tmp)
    }
    names(otu.day1) = paste0("vessel_", 1:num.vessels[i], "_day1")
    names(otu.day14) = paste0("vessel_", 1:num.vessels[i], "_day1")
    
    Y = cbind(otu.day1[,-1],otu.day14[,-1])
    rownames(Y) = seq.names
      
    X.day1 = rbind(diag(num.vessels[i]), rep(0,num.vessels[i]))
    names(X.day1) = paste0("vessel_", 1:num.vessels[i], "_day1")
    X.day14 = rbind(diag(num.vessels[i]), rep(1,num.vessels[i]))
    names(X.day14) = paste0("vessel_", 1:num.vessels[i], "_day14")
    X = cbind(X.day1, X.day14)
    rownames(X) = c(paste0("base::I(vessel)", 1:num.vessels[i]), "base::I(day)14")
    
    upsilon <- nrow(Y) + 3
    Theta <- matrix(0, nrow(Y)-1, nrow(X))
    Gamma <- driver::bdiag(10*diag(num.vessels[i]), 1) # Learn this function as well -- makes block diagonal matricies
    Omega = diag(nrow(Y))
    G = cbind(diag(nrow(Y)-1), -1)
    Xi = G%*%Omega%*%t(G)

    fit_gs <- ssrv.mln(as.matrix(Y),X, covariate = rownames(X)[7],
                                  upsilon = upsilon, Gamma = Gamma,
                                  Theta = Theta, Omega = Omega, Xi = Xi, n_samples = 2000, total_model = "flow", sample.totals =flow_filtered_rep, sample = c(paste0( 1:num.vessels[i], "_day1"),paste0(1:num.vessels[i], "_day14")))
    
    
    fit_unif <- ssrv.mln(as.matrix(Y),X, covariate = rownames(X)[7],
                                    upsilon = upsilon, Gamma = Gamma,
                                    Theta = Theta, Omega = Omega, Xi = Xi, n_samples = 2000,  total_model = "logNormal",sd_lnorm = sqrt(1^2/2), mean_lnorm= rep(c(log(1),log(4)), each = num.vessels[i]),  sample = 1:(2*num.vessels[i]))
    
    sig_tram_gs = sig_tram(fit_gs)
    truth.pos = sig_tram_gs$category
    
    sig_tram_unif = sig_tram(fit_unif)
    
    tmp.fdr[1] = tmp.fdr[1] + ifelse(nrow(sig_tram_unif) > 0, sum(!(sig_tram_unif$category %in% truth.pos))/nrow(sig_tram_unif), 0)
    
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
      rownames_to_column("category") %>% 
      dplyr::select(category, log2FoldChange, padj, lfcSE) %>% 
      mutate(low = log2FoldChange -1.96*lfcSE, 
             high = log2FoldChange + 1.96*lfcSE) %>% 
      mutate(mean=log2FoldChange) %>%
      filter(padj < .05)
    tmp.fdr[3] = tmp.fdr[3] + ifelse(nrow(deseq_results) > 0, sum(!(deseq_results$category %in% truth.pos))/nrow(deseq_results), 0)
    
    ###Aldex2 now
    
    mmX <- t(X)
    aldex_fit <- aldex(Y,mmX,"glm",mc.samples = 128, denom = "all")
    
    aldex_results = aldex_fit %>% 
      rownames_to_column("category") %>%
      dplyr::select(category,model.base..I.day.14.t.value, model.base..I.day.14.Pr...t.., model.base..I.day.14.Pr...t...BH) %>%
      mutate(pval = model.base..I.day.14.Pr...t..) %>%
      mutate(padj = model.base..I.day.14.Pr...t...BH) %>%
      filter(pval < .05)
    tmp.fdr[2] = tmp.fdr[2] + ifelse(nrow(aldex_results) > 0, sum(!(aldex_results$category %in% truth.pos))/nrow(aldex_results), 0)
  }
  fdr[i,2] = tmp.fdr[1]/replicates
  fdr[i,3] = tmp.fdr[2]/replicates
  fdr[i,4] = tmp.fdr[3]/replicates
  
  print(i)
}


fdr.all = data.frame(vals = rep(num.vessels,3), fdr = c(c(fdr[,3]),c(fdr[,4]), c(fdr[,2])), method = rep(c("ALDEx2", "DESeq2", "SSRV (Design)"), each = length(num.vessels)))

fdr.all$method = as.factor(fdr.all$method)

ggplot(fdr.all, aes(x=vals, y=fdr, color=method, fill = method, linetype = method)) +
  geom_line(alpha = 1, lwd = 1.1) +
  theme_bw()+
  xlim(c(6,50)) +
  ylim(c(0,1))+
  xlab("Number of Vessels per Condition") +
  ylab("False Discovery Rate") + 
  scale_color_manual(values=c("#AF4BCE", "#EA7369", "#2b83ba")) + 
  scale_linetype_manual(values = c("dotted", "twodash", "longdash")) +
  theme(text=element_text(size=14)) + 
  geom_hline(yintercept=15/32, linetype="dashed", color = "grey") +
  theme(legend.title = element_blank()) +
  theme(legend.position = c(.65, .825))
ggsave(file.path("results", "unacknowledged_bias_realData.pdf"), height=4, width=4.5)

