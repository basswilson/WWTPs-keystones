##################################################
## Project: WWTPs keystones
## File name: co-occurrence network.R
## Date: 2024-11-04
## Author: Wilson
## R_Version: R version 4.4.1 (2024-06-14 ucrt)
##################################################

library(WGCNA)
library(psych)
library(reshape2)
library(igraph)

#读取相关性网络涵盖的所有OTU
otu = read.csv('otu_table_filt.csv', header = TRUE,row.names = 1)
otu = as.data.frame(t(otu))
#读取metadata
group = read.csv("metadata.csv", header = TRUE)
group[group == "NA"] <- NA
colnames(group)[1] <- "samples"

otu = cbind(rownames(otu),otu)
colnames(otu)[1] = "samples"

df <- merge(otu,group,by="samples")

#以分组“温度10~15”为例
DO5_10 <- otu[which(df$Tem == "10~15"), ]
rownames(DO5_10) <- DO5_10[,1]
DO5_10 <- DO5_10[,-1]

# 计算每一列的和
sums <- colSums(DO5_10)

# 找到和为零的列,判断是否需要删除
a <- which(sums == "0")
a
DO5_10 <- DO5_10[,-a]





cor = corAndPvalue(DO5_10,y=NULL,use = "pairwise.complete.obs", 
                   alternative='two.sided',method='spearman')  #OTU之间的Spearman相关系数和p值



r = cor$cor # 获取相关系数
p = cor$p                       #获取p值
p = p.adjust(p, method = 'BH')  #对p值进行BH校正

r[p > 0.05 | abs(r) < 0.6] = 0  # 对相关性进行筛选，删除假阳性结果：p值>0.001或|r|<0.4的将被去除（赋0值）

write.csv(data.frame(r, check.names = FALSE), 'corr.matrix.csv')     # 将相关系数矩阵写入csv文件中          

g = graph_from_adjacency_matrix(r,mode="undirected",weighted=TRUE,diag = FALSE) #根据相关系数矩阵创建一个加权无向图
g
g = delete.vertices(g, names(degree(g)[degree(g) == 0]))        #删除度数为0的孤立节点

E(g)$corr = E(g)$weight            #为网络的边属性赋值（权重）
E(g)$weight = abs(E(g)$weight)     #为网络的边属性赋值（权重）

#读取taxonomy信息
tax = read.csv('otu注释.csv', row.names=1, header=T)     #读取节点分类信息  

tax = tax[as.character(V(g)$name), ]                  #为节点加上分类信息
V(g)$Kingdom = tax$Kingdom                            #界
V(g)$Phylum = tax$Phylum                              #门
V(g)$Class = tax$Class                                #纲
V(g)$Order = tax$Order                                #目
V(g)$Family = tax$Family                              #科
V(g)$Genus = tax$Genus                                #属
# V(g)$Species = tax$Species                            #种

node_list = data.frame(
  label = names(V(g)),
  kingdom = V(g)$Kingdom,
  phylum = V(g)$Phylum,
  class = V(g)$Class,
  order = V(g)$Order,
  family = V(g)$Family,
  genus=V(g)$Genus)
# species = V(g)$Species)                              #创建节点列表

head(node_list)
write.csv(node_list, 'network.node_list.csv')          #并将其写入csv文件中

edge = data.frame(as_edgelist(g))                      #创建边列表
edge_list = data.frame(
  source = edge[[1]],
  target = edge[[2]],
  weight = E(g)$weight,
  correlation = E(g)$corr
)

head(edge_list)
write.csv(edge_list, 'network.edge_list.csv')          #并将其写入csv文件中

write.graph(g, 'network2.graphml', format = 'graphml')  #后续在Gephi中可视化                         


######计算网络常用的几种拓扑系数#####
nodes_num = length(V(g))                   #节点数
nodes_num

edges_num = length(E(g))                   #边数
edges_num

positive.cor_num = sum(E(g)$corr>0)        #正相关的数量
positive.cor_num

negative.cor_num = sum(E(g)$corr<0)        #负相关的数量
negative.cor_num

average_degree = mean(degree(g))           #平均度
average_degree

average_path_length = average.path.length(g, directed = FALSE)     #平均路径长度
average_path_length

network_diameter = diameter(g, directed = FALSE)                   #网络直径
network_diameter

network_density = graph.density(g)                                 #网络密度
network_density

clustering_coefficient = transitivity(g)                           #聚类系数
clustering_coefficient


community_result = cluster_fast_greedy(g, 
                                       modularity = TRUE, 
                                       membership = TRUE)          
membership_vector = community_result$membership
modularity_score = modularity(g, 
                              membership = membership_vector)      #计算模块性
modularity_score


network_parameter = data.frame(nodes_num, 
                               edges_num, 
                               positive.cor_num, 
                               negative.cor_num, 
                               average_degree,
                               average_path_length,
                               network_diameter, 
                               network_density,
                               clustering_coefficient,
                               modularity_score
)

network_parameter
write.csv(network_parameter, 'network_parameter.csv')                                  

otu1 = otu
otu1[otu1>0] = 1

write.csv(otu1, 'adjacent_matrix.csv')

