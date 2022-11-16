####Sampling functions##############################################################
create_true_abundances <- function(d, n){
  dd <- length(d)/2
  dat <- d %>%
    sapply(function(x) rpois(n, lambda=x)) %>% 
    t() %>%
    as.data.frame() %>%
    split(rep(1:2, each=dd)) %>%
    purrr::map(~`rownames<-`(.x, paste0("Taxa", 1:dd))) %>%
    purrr::map(t) %>%
    do.call(rbind, .) %>%
    as.data.frame() %>%
    cbind(Condition=factor(rep(c("Pre", "Post"), each=n), levels = c("Pre", "Post")), .) %>%
    `rownames<-`(., NULL)
  return(dat)
}

resample_data <- function(dat, seq.depth){
  ddat <- driver::miniclo(as.matrix(dat[,-1]))
  for (i in 1:nrow(dat)){
    dat[i,-1] <- rmultinom(1, size=seq.depth, prob=ddat[i,])
  }
  return(dat)
}

####Utility functions###############################################################
##Geometric mean
gm <- function(x, na.rm = TRUE){
  exp(mean(log(x[x > 0]), na.rm=na.rm))
}

append_sig <- function(fit, f){
  s <- f(fit)$category
  mutate(fit, sig=ifelse(category%in%s, TRUE, FALSE))
}

strip_end_num <- function(x){
  x <- as.character(x)
  as.numeric(str_extract(x, "[[:digit:]]+"))
}

sig_code <- function(sig, Taxa, truth){
  out <- rep("TN", length(Taxa))
  out[sig &(Taxa %in% truth)] <- "TP" # True Positives
  out[sig & (out!="TP")] <- "FP" # False Positives
  out[!sig & (Taxa %in% truth)] <- "FN" # False Negatives
  return(out)
}

####Plotting functions##############################################################
plot_count <- function(dat){
  gather(dat, Taxa, Count, -Condition) %>% 
    mutate(Taxa=strip_end_num(Taxa)) %>% 
    mutate(Taxa=factor(Taxa)) %>% 
    ggplot(aes(x=Taxa, y=Count)) +
    geom_boxplot(aes(fill = Condition, color = Condition), position=position_dodge(width=1), 
                size=1)+
    scale_y_log10() +
    theme(text = element_text(size=16)) +
    theme_bw() +
    scale_fill_manual(values = c("#fdae61", "#2b83ba")) + 
    scale_color_manual(values = c("#fdae61", "#2b83ba")) +
    labs(color='Antibiotic\nTreatment') +
    labs(fill='Antibiotic\nTreatment') 
}

plot_sig2 <- function(rrs, truth, ...){
  rrs$dat <- NULL
  names(rrs) <- model.names[names(rrs)]
  graph_df = bind_rows(rrs, .id="Model") %>% 
    dplyr::select(Model, category, sig) %>% 
    mutate(Taxa = category) %>% 
    mutate(Taxa=strip_end_num(Taxa)) %>% 
    mutate(sigcode = sig_code(sig, Taxa, truth)) %>% 
    mutate(Taxa=factor(Taxa), sigcode=factor(sigcode, 
                                             levels=c("TP", "TN", 
                                                      "FP", "FN"))) %>% 
    mutate(Model=factor(Model, levels=model.name.levels))
  rownames(graph_df) <- 1:nrow(graph_df)
  ggplot(graph_df, aes(x=Taxa, y=Model)) +
    geom_tile_pattern(aes(fill=sigcode, pattern = sigcode), color="darkgrey",pattern_fill = 'grey',pattern_colour  = 'grey', pattern_density = 0.015) +
    theme_minimal() +
    theme(panel.grid = element_blank(), 
          legend.title=element_blank(),
          text = element_text(size=16),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    scale_pattern_manual(values = c(TP = "none", TN = "none", FP = "none", FN = "stripe")) +
    scale_fill_manual(values= c("black", "white", "grey", "white"))
}

