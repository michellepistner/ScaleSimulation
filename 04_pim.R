# PIM prior
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
                    Theta = Theta, Omega = Omega, Xi = Xi, Theta.t = Theta.t, n_samples = 1000, total_model = "pim", prob = 0.05, flag = TRUE)
##Extracting the total
totals <- apply(exp(fit_pim), MARGIN = c(3), colSums)
theta.perp <- apply(log(totals), 2, FUN= function(vec) mean(vec[7:12]) - mean(vec[1:6]))
