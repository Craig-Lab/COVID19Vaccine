library(ggpubr)
library(ggplot2)
library(R.matlab)
library(corrplot)
library(dplyr)
library(readxl)
setwd("/Users/dengxiaoyan/Desktop/PhD projects/Humoral immune response/Paper_new/Code_new/Matlab_new")
hcw_par <- read_xlsx("hcw_neu_par.xlsx", sheet = 1)#sheet1 combined1; sheet2 combined2; sheet3 constant
senior_par <- read_xlsx("senior_neu_par.xlsx", sheet = 1)
hcw_par$Group <- rep('A',length(hcw_par$id))
senior_par$Group <- rep('B',length(senior_par$id))
all_par <-rbind(hcw_par,senior_par)

wilcox.test(E0 ~ Group, data = all_par, exact = FALSE)
wilcox.test(EC50 ~ Group, data = all_par, exact = FALSE)
wilcox.test(gamma ~ Group, data = all_par, exact = FALSE)

cor.test(all_par$Age, all_par$E0, method = c("pearson", "kendall", "spearman"))
cor.test(all_par$Age, all_par$EC50, method = c("pearson", "kendall", "spearman"))
cor.test(all_par$Age, all_par$gamma, method = c("pearson", "kendall", "spearman"))



p2 <- ggboxplot(
  all_par, x = "Group", y = "E0",
  width = 0.6, color = "Group", fill = "Group",
  palette = c("#b24a4b", "#4276bc"), alpha = 0.2,
  add = "jitter", add.params = list(size = 0.3)
) +
  scale_x_discrete(
    labels = c('HCW', 'Senior')  # Changed 'label' to 'labels'
  )  +
  theme(
    plot.margin = unit(c(0, 0, 0, 0), "inches"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "none",
    axis.line = element_line(linewidth = 0.1),
    axis.text = element_text(size = 8)
  ) +
  scale_y_continuous(
    breaks = c( 0.2,0.3, 0.4, 0.5),  # Specify the desired tick marks
    limits = c(0.1, 0.5)      # Ensure the limits are consistent
  )#+ scale_y_continuous(limits = c(0,1), breaks = c(0.2,0.5,0.8))#+ stat_compare_means(method = "t.test", aes(label = sprintf("p = %5.4f", as.numeric(..p.format..))), size = 5, label.y = 0.8, label.x = 1.6)
p2$layers[[1]]$aes_params$linewidth <- 0.3  # Set this to a finer value if needed
ggsave("boxplot_E0.pdf", width = 2.8 , height = 3.6, units = "cm")

p3 <- ggboxplot(
  all_par, x = "Group", y = "EC50" , width = 0.6, color = "Group", fill= "Group" , palette = c("#b24a4b", "#4276bc"), alpha=0.2, add="jitter",
  add.params = list(size = 0.3) ) +
  scale_x_discrete(
    label = c('HCW','Senior')
  )+
  scale_y_continuous(
    breaks = c(1, 3, 5),  # Specify the desired tick marks
    limits = c(1, 5)      # Ensure the limits are consistent
  ) +
  theme(
    plot.margin = unit(c(0, 0, 0, 0), "inches"),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(), 
    #axis.ticks.y = c('0','0.5','1'),
    legend.position="none",
    axis.line = element_line(linewidth = 0.1),
    axis.text=element_text(size=8)
  ) #+ scale_y_continuous(limits = c(0,7))#+ stat_compare_means(method = "t.test", aes(label = sprintf("p = %5.4f", as.numeric(..p.format..))), size = 5, label.y = 9, label.x = 1.6)
p3
p3$layers[[1]]$aes_params$linewidth <- 0.3  # Set this to a finer value if needed
#compare_means(EC50 ~ Group, data = all_par, method = "wilcox.test")
ggsave("boxplot_EC50.pdf",  width = 2.8 , height = 3.6,   units = "cm")

p1 <- ggboxplot(
  all_par, x = "Group", y = "gamma" , , width = 0.6, color = "Group", fill= "Group" , palette = c("#b24a4b", "#4276bc"), alpha=0.2, add="jitter",
  add.params = list(size = 0.3) )+  scale_x_discrete(
  label = c('HCW','Senior')
) +
  theme(
    plot.margin = unit(c(0, 0, 0, 0), "inches"),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(), 
    legend.position="none",
    axis.line = element_line(linewidth = 0.1),
    axis.text=element_text(size=8)
  )+
  scale_y_continuous(
    breaks = c(2,3, 4),  # Specify the desired tick marks
    limits = c(2, 4)       # Ensure the limits are consistent
  )  #+ stat_compare_means(method = "t.test", aes(label = sprintf("p = %5.4f", as.numeric(..p.format..))), size = 5, label.y = 5.5, label.x = 1.6)
p1$layers[[1]]$aes_params$linewidth <- 0.3  # Set this to a finer value if needed
p1
ggsave("boxplot_h.pdf",   width = 2.8 , height = 3.6,   units = "cm")




## plot correlation matrix
library(ggplot2)
library(readxl)
library(reshape2)

M <- cor(all_par[, c("Age", "E0", "EC50", "gamma")])


# Custom function to calculate p-values for correlation matrix
cor.mtest <- function(mat, conf.level = 0.95) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], conf.level = conf.level)
      p.mat[i, j] <- tmp$p.value
      p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  return(p.mat)
}

# Calculate p-values
p_values <- cor.mtest(all_par[, c("Age", "E0", "EC50", "gamma")])

corr_data <- melt(M)
pval_data <- melt(p_values)
colnames(corr_data) <- c("Var1", "Var2", "Correlation")
colnames(pval_data) <- c("Var1", "Var2", "p_value")


plot_data <- merge(corr_data, pval_data)

plot_data$Var2 <- factor(plot_data$Var2, levels = c("gamma", "EC50", "E0", "Age"))

plot_data$stars <- cut(plot_data$p_value,
                       #breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                       breaks = c(-Inf, 0.05, Inf),
                       #labels = c("***", "**", "*", ""),
                       labels = c( "*", ""),
                       right = FALSE)
col_palette <- colorRampPalette(c('#4276bc', '#7896c1', '#9baecb', '#bdc6d7', '#faf5f4', '#e1b9b6', '#d39794', '#c67875', '#b24a4b'))

# Create and save the plot
pdf(file = "correlation_neu_age_plot.pdf", width = 4 / 2.54, height = 4 / 2.54)
ggplot(data = plot_data, aes(x = Var1, y = Var2, fill = Correlation)) +
  geom_tile(color = "grey") +  # Add grey border for each cell
  scale_fill_gradientn(colors = col_palette(10), limits = c(-1, 1)) +
  geom_text(aes(label = round(Correlation, 2)), color = "black", size = 8 / .pt, hjust = 0.5, vjust = 0.5) +
  geom_text(aes(label = stars), color = "black", size = 7 / .pt, hjust = 0.5, vjust = -0.5) +  # Stars above numbers
  theme_minimal() +
  theme(
    axis.text = element_blank(),       # Remove axis text
    axis.title = element_blank(),      # Remove axis titles
    axis.ticks = element_blank(),      # Remove axis ticks
    panel.grid = element_blank(),
    panel.border = element_blank(),
    legend.position = "none"  ) +
  coord_fixed()  # Ensure the plot stays square
dev.off()




