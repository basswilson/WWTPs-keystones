##################################################
## Project: WWTPs keystones
## File name: α diversity.R
## Date: 2024-11-04
## Author: Wilson
## R_Version: R version 4.4.1 (2024-06-14 ucrt)
##################################################


library(vegan)
library(ggplot2)
library(tidyverse)
library(agricolae)
library(car)
library(FSA)


#读取数据，保存在dune和dune.env中
dune <- read.csv("61个基石物种的丰度.csv", row.names = 1)
dune.env <- read.csv("metadata.csv", row.names = 1)
a<-which(!is.na(dune.env$DO))
dune <- dune[a,]
dune.env <- dune.env[a,]

#计算alpha多样性，包括均匀度、丰富度和香农指数。

#计算丰富度richness，即群落中丰度大于0的otu数量之和
richness <- rowSums(dune>0)
#计算香农指数Shannon diversity index，以e作为底数
shannon <- diversity(dune, index = 'shannon', base = exp(1))
#计算均匀度Pielou evenness，即香农指数与ln(Richness)的比值
pielo <- shannon/log(richness, base = exp(1))

#创建函数一步计算alpha多样性。
calculate_alpha <- function(otu){
  data_richness <- rowSums(otu>0)
  data_shannon <- diversity(otu, index = 'shannon', base = exp(1))
  data_pielou <- data_shannon/log(data_richness, base = exp(1))
  alpha_matrix <- cbind(data_pielou, data_richness, data_shannon)
  alpha_df <- as.data.frame(alpha_matrix)
  return(alpha_df)
}

#利用函数calculate_alpha()，创建绘图所需数据框并重命名其中的列：
plot_df <- calculate_alpha(dune)%>%
  cbind(dune.env$DO)%>%
  rename_with(~"Group", 4)


#绘制Shannon的箱线图
p <- p_shannon <- ggplot(plot_df, aes(Group, data_shannon))+
  geom_boxplot(aes(fill = Group))

p

# ###############分界线##################分界线##################分界线###########
# 
# 
# 
# #到这一步，基本的分面箱线图已经完成，但我们注意到原图做了差异分析，并用字母标记了差异分析结果。
# #首先验证总体是否符合正态分布。
# 
# #执行正态性检验（Shapiro-Wilk检验）
# # 假设数据存储在变量data_shannon中
# # H0(原假设): 数据符合正态分布
# # Ha(备选假设): 数据不符合正态分布
# shapiro_test <- shapiro.test(plot_df$data_shannon)
# print(shapiro_test)
# 
# #如果P小于显著性水平，即p<0.5，我们可以拒绝原假设，
# #并得出结论：在给定的显著性水平下，数据_shannon在该样本中不服从正态分布。
# 
# #执行方差齐性检验(Levene检验)
# # 假设数据存储在变量data_shannon中，按照组别存储在变量Group中
# # H0: 不同组的数据方差相等（方差齐性）
# # Ha: 不同组的数据方差不等
# levene_test <- leveneTest(data_shannon ~ Group, data = plot_df)
# levene_test
# 
# #如果如果P小于显著性水平，即p<0.5，我们可以拒绝原假设，
# #并得出结论：在给定的显著性水平下，数据_shannon在该样本中不服从方差齐性。
# 
# 
# 
# 
# 
# 
# 
# 
# #根据以上结果有两种方法来分析，跳转到对应方法的代码继续跑即可。
# 
# ##################方法一#################方法一################################
# 
# #方法一 若满足正态分布和方差齐性，则进行one-way ANOVA分析，代码如下。
# #先选择Shannon进行one-way ANOVA分析
# shannon_anova <- aov(data_shannon~Group, data = plot_df)
# #查看ANOVA结果，其中Pr(>F)是p值，p<0.05时，认为不同组的Shannon具有显著差异
# summary(shannon_anova)
# #进一步得到更详细的两两比较，用TukeyHSD()函数
# pair_comparison <- TukeyHSD(shannon_anova)
# pair_comparison <- as.data.frame(pair_comparison$Group)
# pair_comparison
# 
# #计算分组平均数
# group_mean <- aggregate(x = plot_df$data_shannon, by = list(plot_df$Group), FUN = mean)%>%
#   rename_with(~c("Group", "mean_val"), 1:2)
# 
# 
# #创建一个pvalue矩阵
# ntr <- nrow(group_mean)
# mat <- matrix(1, ncol = ntr, nrow = ntr)
# p <- pair_comparison$`p adj`
# k <- 0
# for (i in 1:(ntr - 1)) {
#   for (j in (i + 1):ntr) {
#     k <- k + 1
#     mat[i, j] <- p[k]
#     mat[j, i] <- p[k]
#   }
# }
# 
# 
# treatments <- as.vector(group_mean$Group)
# means <- as.vector(group_mean$mean_val)
# alpha <- 0.05
# pvalue <- mat
# out <- orderPvalue(treatments, means, alpha, pvalue, console = TRUE)
# out
# 
# 
# #加入字母差异，做图。
# p <- p_shannon <- ggplot(plot_df, aes(Group, data_shannon))+
#   geom_boxplot(aes(fill = Group))+
#   geom_text(data = out, aes(x = rownames(out), y = 0, label = groups),
#             vjust = 1, color = "black")+  # 调整 x 和 y 的位置和 vjust 的值以使文本正确显示
#   theme(panel.grid = element_blank(), 
#         panel.background = element_rect(fill = 'white'), 
#         panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid"), 
#         axis.title = element_blank(), 
#         axis.text.x = element_blank())
# 
# #输出，文件名为p_shannon_ANOVA.svg
# ggsave("p_shannon.svg", p, width = 8, height = 8)




##################方法二#################方法二################################

#方法二若不满足正态分布和方差齐性，特别是不满足正态分布，则进行Kruskal-Wallis检验

#如果p值小于预设的显著性水平（通常为0.05）
#则可以拒绝原假设，认为不同组别之间存在显著差异。
#进一步得到更详细的两两比较，用dunn()函数
# 执行Kruskal-Wallis检验
kw_test <- kruskal.test(data_shannon ~ Group, data = plot_df)
print(kw_test)
pvalue <- kw_test$p.value
pvalue <- sprintf("%.2e", kw_test$p.value)



#加入字母差异，做图。
p <- p_shannon <- ggplot(plot_df, aes(Group, data_shannon))+
  geom_boxplot(aes(fill = Group))+
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(fill = 'white'), 
        panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid"), 
        axis.title = element_blank(), 
        axis.text.x = element_blank())
p <- p +
  annotate("text", x = Inf, y = Inf, 
           label = paste("Kruskal-Wallis test, p =", pvalue), 
           hjust = 1, vjust = 4, size = 3, parse = FALSE)
p

ggsave("DO.svg", plot = p, width = 4, height = 4)
