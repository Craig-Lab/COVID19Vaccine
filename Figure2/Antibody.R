library(ggpubr)
library(ggplot2)
library(R.matlab)
library(corrplot)
library(dplyr)
library(readxl)
setwd("/Users/dengxiaoyan/Desktop/PhD projects/Humoral immune response/Paper_new/Code_new/Matlab_new")
hcw_par <- read_xlsx("hcw_ant_par.xlsx", sheet = 1)
senior_par <- read_xlsx("senior_ant_par.xlsx", sheet = 1)
hcw_par$Group <- rep('HCW',length(hcw_par$id))
senior_par$Group <- rep('Senior',length(senior_par$id))
all_par <-rbind(hcw_par,senior_par)
wilcox.test(d_t ~ Group, data = all_par, exact = FALSE)
wilcox.test(lam2 ~ Group, data = all_par, exact = FALSE)
wilcox.test(lam3 ~ Group, data = all_par, exact = FALSE)

wilcox.test(1+lam2 ~ Group, data = all_par, exact = FALSE)
wilcox.test(1+lam2+lam3 ~ Group, data = all_par, exact = FALSE)

cor.test(all_par$Age, all_par$d_t, method = c("pearson", "kendall", "spearman"))
cor.test(all_par$Age, all_par$lam2, method = c("pearson", "kendall", "spearman"))
cor.test(all_par$Age, all_par$lam3, method = c("pearson", "kendall", "spearman"))
summary_stats <- all_par %>%
  group_by(Group) %>%
  summarise(across(where(is.numeric), 
                   list(Min = min, 
                        Max = max, 
                        Median = median, 
                        SD = sd), 
                   na.rm = TRUE))

# Display the results
print(summary_stats)



p4 <- ggboxplot(
  all_par, x = "Group", y = "d_t" ,  width = 0.6, color = "Group", fill= "Group" , palette = c("#b24a4b", "#4276bc"), alpha=0.2, add="jitter",
  add.params = list(size = 0.3) ) + scale_x_discrete(
  label = c('HCW','Senior')
) +
  theme(
    plot.margin = unit(c(0, 0, 0, 0), "inches"),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(), 
    legend.position="none",
    axis.line = element_line(linewidth = 0.1),
    axis.text=element_text(size= 8)
  ) + 
  scale_y_continuous(
    breaks = c( 0.5,1.0),  # Specify the desired tick marks
    limits = c(0.25, 1.25)      # Ensure the limits are consistent
  )
# Reduce boxplot line thickness
p4$layers[[1]]$aes_params$linewidth <- 0.3  # Set this to a finer value if needed
p4
ggsave("boxplot_dt.pdf", width = 3.2 , height = 3.6, units = "cm")


p5 <- ggboxplot(
  all_par, x = "Group", y = "lam2" , width = 0.6, color = "Group", fill= "Group" , 
  palette = c("#b24a4b", "#4276bc"), alpha = 0.2, add = "jitter",
  add.params = list(size = 0.3)) +
  scale_x_discrete(
    labels = c('HCW', 'Senior')
  ) +
  scale_y_continuous(
    breaks = c(5, 10, 15),  # Specify the desired tick marks
    limits = c(5, 15)       # Ensure the limits are consistent
  ) +
  theme(
    plot.margin = unit(c(0, 0, 0, 0), "inches"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "none",
    axis.line = element_line(linewidth = 0.1),
    axis.text = element_text(size = 8)
  )

p5$layers[[1]]$aes_params$linewidth <- 0.3  # Set this to a finer value if needed

p5
ggsave("boxplot_lam2.pdf", width = 3.2 , height = 3.6, units = "cm")



p6 <- ggboxplot(
  all_par, x = "Group", y = "lam3" , width =  0.6,  color = "Group", fill= "Group" , palette = c("#b24a4b", "#4276bc"), alpha=0.2, add="jitter",
  add.params = list(size = 0.3))+ scale_x_discrete(
  label = c('HCW','Senior')
) +
  theme(
    plot.margin = unit(c(0, 0, 0, 0), "inches"),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    legend.position="none",
    axis.line = element_line(linewidth = 0.1),
    axis.text=element_text(size = 8)
  )+ ylim(4.5, 20)
p6$layers[[1]]$aes_params$linewidth <- 0.3  # Set this to a finer value if needed

p6
ggsave("boxplot_lam3.pdf", width = 3.2, height = 3.6, units = "cm")


# Load necessary libraries
library(ggplot2)
library(reshape2)
# Calculate the correlation matrix
M <- cor(all_par[, c("Age", "d_t", "lam2", "lam3")])

# Custom function to calculate p-values for correlation matrix
cor.mtest <- function(mat, conf.level = 0.95) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], conf.level = conf.level)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  return(p.mat)
}

# Generate the p-value matrix
p_values <- cor.mtest(all_par[, c("Age", "d_t", "lam2", "lam3")])
corr_data <- melt(M)
pval_data <- melt(p_values)
colnames(corr_data) <- c("Var1", "Var2", "Correlation")
colnames(pval_data) <- c("Var1", "Var2", "p_value")


plot_data <- merge(corr_data, pval_data)

plot_data$Var2 <- factor(plot_data$Var2, levels = c("lam3", "lam2", "d_t", "Age"))

# Define significance stars based on p-values
plot_data$stars <- cut(plot_data$p_value,
                       breaks = c(-Inf,  0.05, Inf),
                       #labels = c("***", "**", "*", ""),
                       labels = c("*", ""),
                       right = FALSE)

# Custom color palette for correlation values
col_palette <- colorRampPalette(c('#4276bc', '#7896c1', '#9baecb', '#bdc6d7', '#faf5f4', '#e1b9b6', '#d39794', '#c67875', '#b24a4b'))

pdf(file = "correlation_ant_age.pdf", width = 4/ 2.54, height = 4 / 2.54)

ggplot(data = plot_data, aes(x = Var1, y = Var2, fill = Correlation)) +
  geom_tile(color = "grey") +  # Add grey border for each cell
  scale_fill_gradientn(colors = col_palette(10), limits = c(-1, 1)) +
  geom_text(aes(label = round(Correlation, 2)), color = "black", size = 8 / .pt, hjust = 0.5, vjust = 0.5) +
  geom_text(aes(label = stars), color = "black", size = 7 / .pt, hjust = 0.5, vjust = -0.5) +  # Stars above numbers
  theme_minimal() +
  theme(
    axis.text = element_blank(),       # Remove axis text (parameter names)
    axis.title = element_blank(),      # Remove axis titles
    axis.ticks = element_blank(),      # Remove axis ticks
    panel.grid = element_blank(),      # Remove gridlines
    panel.border = element_blank(),
       legend.position = "none"           # Remove legend
  )
dev.off()
# Plot with customized appearance



