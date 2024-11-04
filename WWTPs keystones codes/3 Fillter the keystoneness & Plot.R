##################################################
## Project: WWTPs keystones
## File name: Fillter the keystoneness & Plot.R
## Date: 2024-11-04
## Author: Wilson
## R_Version: 4.4.1 (2024-06-14 ucrt)
##################################################

# 读取keystoness表格
keystoness <- read.csv("keystoness.csv", header = TRUE)

# 读取species_ID跟属名对应的表
id <- read.csv("ID跟属对应.csv",header = T)


# 根据每个species的str_pred计算中位数。
# 加载包
library(dplyr)

# 计算每个 species_id 对应的 str_pred 中位数
median_str_pred <- keystoness %>%
  group_by(species_id) %>%
  summarize(median_str_pred = median(str_pred))


# 计算每个 species_id 对应的 str_true中位数
median_str_true <- keystoness %>%
  group_by(species_id) %>%
  summarize(median_str_true = median(str_true))

# 合并计算出来的中位数
median <- merge(median_str_pred, median_str_true, by = "species_id", all.x = TRUE)

################# 随后根据Top前20绘制山峦图    ################

# 按照名为median_str_pred的列重新对整个数据框降序
median_str_pred <- median_str_pred %>%
  arrange(desc(median_str_pred))

# 筛选Top 20
median_str_pred <- median_str_pred[-1,]
top_20 <- median_str_pred[1:20,]
# 按照ID跟属对应，合并表格到top_20中
top_20 <- merge(top_20, id, by = "species_id", all.x = TRUE)


top_20 <- top_20 %>%
  arrange(desc(median_str_pred))

# 过滤
# 名为 top_20 的数据框，包含了你想要保留的前 20个species_id
# 确保 top_20 数据框中有一个名为 species_id 的列

# 根据Top20中位数，将所有的species的keystoness筛选，方便之后的山峦图作图
# 使用 filter() 函数根据 top_20 数据框的 species_id 列对 keystoness 数据框进行过滤
filtered_keystoness_top_20 <- filter(keystoness, species_id %in% top_20$species_id)
filtered_keystoness_top_20 <- merge(filtered_keystoness_top_20, id, by = "species_id", all.x = TRUE)
filtered_keystoness_top_20$genus <- factor(filtered_keystoness_top_20$genus, levels = unique(top_20$genus))
library(dplyr)

filtered_keystoness_top_20 <- filtered_keystoness_top_20 %>% arrange(row_number())

# 绘制山峦图


library(ggridges)
library(ggplot2)
library(cols4all)


# 假设数据框名为keystoness


ggplot(data = filtered_keystoness_top_20,
             aes(x = str_pred, y = genus, fill = genus)) +
  geom_density_ridges(alpha = 0.8,
                      color = 'white',
                      rel_min_height = 0.01,
                      scale = 1.8,
                      quantile_lines = TRUE,
                      quantiles = 2
  ) +
  
  theme_linedraw() +
  theme(legend.position = 'none',
        axis.text.x = element_text(size = 13),  # 设置 x 轴字体大小为 10
        axis.text.y = element_text(size = 13)) +  # 设置 y 轴字体大小为 10
  scale_x_continuous(limits = c(0, 0.05),
                     breaks = seq(0, 0.05, by = 0.01)) +
  scale_y_discrete(limits = rev(levels(filtered_keystoness_top_20$genus)), position = "right") # 反转Y轴顺序，把Y轴位置放在右边


write.csv(top_20, file = "top_20.csv")
ggsave("Top_20.svg", width = 10, height = 6, dpi = 300, limitsize = FALSE)
ggsave("Top_20.png", width = 10, height = 6, dpi = 300, limitsize = FALSE)






################# 随后根据bottom20绘制山峦图 ################

# 按照名为median_str_pred的列重新对整个数据框升序
median_str_pred <- median_str_pred %>%
  arrange(median_str_pred)

# 筛选bottom20
bottom_20 <- median_str_pred[1:20,]
# 按照ID跟属对应，合并表格到top_20中
bottom_20 <- merge(bottom_20, id, by = "species_id", all.x = TRUE)


bottom_20 <- bottom_20 %>%
  arrange(median_str_pred)

# 过滤
# 名为 bottom_20 的数据框，包含了你想要保留的前 20个species_id
# 确保 bottom_20 数据框中有一个名为 species_id 的列

# 根据bottom_20中位数，将所有的species的keystoness筛选，方便之后的山峦图作图
# 使用 filter() 函数根据 bottom_20 数据框的 species_id 列对 keystoness 数据框进行过滤
filtered_keystoness_bottom_20 <- filter(keystoness, species_id %in% bottom_20$species_id)
filtered_keystoness_bottom_20 <- merge(filtered_keystoness_bottom_20, id, by = "species_id", all.x = TRUE)
filtered_keystoness_bottom_20$genus <- factor(filtered_keystoness_bottom_20$genus, levels = unique(bottom_20$genus))
# 绘制山峦图

# 假设数据框名为keystoness


ggplot(data = filtered_keystoness_bottom_20,
       aes(x = str_pred, y = genus, fill = genus)) +
  geom_density_ridges(alpha = 0.8,
                      color = 'white',
                      rel_min_height = 0.01,
                      scale = 1.8,
                      quantile_lines = TRUE,
                      quantiles = 2
  ) +
  
  theme_linedraw() +
  theme(legend.position = 'none',
        axis.text.x = element_text(size = 13),  # 设置 x 轴字体大小为 10
        axis.text.y = element_text(size = 13)) +  # 设置 y 轴字体大小为 10
  scale_x_continuous(limits = c(0, 0.01),
                     breaks = seq(0, 0.01, by = 0.002)) +
  scale_y_discrete(position = "right")
  
#此处就不需要反转Y顺序了

write.csv(bottom_20, file = "bottom_20.csv")
ggsave("bottom_20.png", width = 7, height = 6, dpi = 300, limitsize = FALSE)
ggsave("bottom_20.svg", width = 7, height = 6, dpi = 300, limitsize = FALSE)



