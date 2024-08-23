
# test ExoGRU

x <- "/BioII/lulab_b/baopengfei/projects/motif-RBP-EDA/output/lulab/11RNAsmallDomain_EVenrich_gn.bed.csv"
EV.exo <- read.csv(x)
EV.exo$type <- "EV"

gn.exo <- list()
for (SEED in 100:109){
  # SEED <- 100
  message(SEED)
  x <- paste0("/BioII/lulab_b/baopengfei/projects/motif-RBP-EDA/output/lulab/shuf/11RNAsmallDomain_EVenrich_gn.bed.shuffle",SEED,".csv")
  tmp <- read.csv(x)
  gn.exo[[SEED]] <- tmp
}
gn.exo <- as.data.frame(do.call(rbind,gn.exo))
gn.exo$type <- "gn_bg"

tx.exo <- list()
for (SEED in 100:109){
  # SEED <- 100
  message(SEED)
  x <- paste0("/BioII/lulab_b/baopengfei/projects/motif-RBP-EDA/output/lulab/shuf/11RNAsmallDomain_EVenrich.bed.shuffle",SEED,"_gn.csv")
  tmp <- read.csv(x)
  tx.exo[[SEED]] <- tmp
}
tx.exo <- as.data.frame(do.call(rbind,tx.exo))
tx.exo$type <- "tx_bg"

df <- as.data.frame(rbind(EV.exo,gn.exo,tx.exo))


# plot bar
library(ggplot2)
library(ggpubr)

df$method <- df$type
df$prob
p11 <- ggplot(data = df, aes(x = method, y = prob, fill = type)) +
  geom_bar(stat = "summary", fun = "mean", position = position_dodge(), width = 0.7) +
  geom_errorbar(stat = "summary", fun.data = mean_se, position = position_dodge(0.7), width = 0.2) +
  # facet_wrap(~site, scales = "free", labeller = label_both) +
  stat_compare_means(aes(group = type), method = "wilcox.test", 
                     comparisons = list(c("EV","gn_bg"),c("EV","tx_bg")),
                     # method.args = list(alternative = "less"),
                     label = "p.signif", 
                     label.x.npc = 0.8, label.y.npc = 0.8, 
                     hide.ns = F, size = 10, paired = FALSE
                     ) +
  scale_fill_manual(values = c('firebrick', 'salmon2','salmon1')) +
  labs(title = "", x = "", y = "EV prob. (ExoGRU)") +
  ylim(c(0, 1.4)) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.x = element_text(size = 20, hjust = 1, vjust = 0.5, angle = 90),
    axis.text.y = element_text(size = 20),
    plot.title = element_text(size = 20),
    legend.position = "none",
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16)
  )

print(p11)


