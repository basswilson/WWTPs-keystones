##################################################
## Project: WWTPs keystones
## File name: Cohesion.R
## Date: 2024-11-04
## Author: Wilson
## R_Version: R version 4.4.1 (2024-06-14 ucrt)
##################################################


## 保留物种的阈值。如共有100个样本，0.1为只保留出现在大于10个样本中的物种。
pers.cutoff <- 0.10
## 迭代次数。推荐>= 200，但是速度会慢。 
iter <- 200

tax.shuffle <- T
use.custom.cors <- F

#读入数据，行为样本，列为物种
b <- read.csv("genus20带属名.csv", header = T, row.names = 1)

##三个函数
#1.统计0的个数
zero <- function(vec){
  num.zero <- length(which(vec == 0))
  return(num.zero)
}

#2.挑出来所有相关性的负值并求平均值。
neg.mean <- function(vector){
  neg.vals <- vector[which(vector < 0)]
  n.mean <- mean(neg.vals)
  if(length(neg.vals) == 0) n.mean <- 0
  return(n.mean)
}

#3.挑出来所有相关性的正值并求平均值。
pos.mean <- function(vector){
  pos.vals <- vector[which(vector > 0)]
  p.mean <- mean(pos.vals)
  if(length(pos.vals) == 0) p.mean <- 0
  return(p.mean)
}

# 读入自己的相关性矩阵
if(use.custom.cors == T) {
  custom.cor.mat <- read.csv("your_path_here.csv", header = T, row.names = 1)
  custom.cor.mat <- as.matrix(custom.cor.mat)
  #Check that correlation matrix and abundance matrix have the same dimension
  print(dim(b)[2] == dim(custom.cor.mat)[2])
}

# 去掉都是0的行和列
c <- as.matrix(b)
c <- c[rowSums(c) > 0, colSums(c) > 0]

# 计算样本的总序列数
rowsums.orig <- rowSums(c)

# 确定物种个数的阈值
zero.cutoff <- ceiling(pers.cutoff * dim(c)[1])

# 根据阈值筛选物种
d <- c[ , apply(c, 2, zero) < (dim(c)[1]-zero.cutoff) ]
# 移除丰度都是0的样本
d <- d[rowSums(d) > 0, ]

#如果自己提供了相关性表，也要根据阈值筛选
if(use.custom.cors == T){
  custom.cor.mat.sub <- custom.cor.mat[apply(c, 2, zero) < (dim(c)[1]-zero.cutoff), apply(c, 2, zero) < (dim(c)[1]-zero.cutoff)]
}

# 计算丰度的百分比
rel.d <- d / rowsums.orig

# 计算实际相关矩阵
cor.mat.true <- cor(rel.d)

#计算零模型
med.tax.cors <- vector()
if(use.custom.cors == F) {
  if(tax.shuffle) {  #第一种零模型计算方法，打乱列
    for(which.taxon in 1:dim(rel.d)[2]){
      perm.cor.vec.mat <- vector()
      
      for(i in 1:iter){ 
        perm.rel.d <- matrix(numeric(0), dim(rel.d)[1], dim(rel.d)[2])
        rownames(perm.rel.d) <- rownames(rel.d)
        colnames(perm.rel.d) <- colnames(rel.d)
        
        #对每个物种
        for(j in 1:dim(rel.d)[2]){  
          # 先替换所有的值
          perm.rel.d[, j ] <- sample(rel.d[ ,j ]) 
        }
        
        # focal物种维持原状不随机
        perm.rel.d[, which.taxon] <- rel.d[ , which.taxon]
        
        # 算相关性
        cor.mat.null <- cor(perm.rel.d)
        
        # 保存每次对于focal物种迭代的结果
        perm.cor.vec.mat <- cbind(perm.cor.vec.mat, cor.mat.null[, which.taxon])
        
      }
      med.tax.cors <- cbind(med.tax.cors, apply(perm.cor.vec.mat, 1, median))
      # 运行情况
      if(which.taxon %% 20 == 0){print(which.taxon)}
    }
  } else { #第二种零模型。打乱样本内所有丰度，打乱行
    for(which.taxon in 1:dim(rel.d)[2]){  
      perm.cor.vec.mat <- vector()
      
      for(i in 1:iter){ 
        perm.rel.d <- rel.d 
        #对每个物种
        for(j in 1:dim(rel.d)[1]){ 
          #先取出大于0的值
          which.replace <- which(rel.d[j, ] > 0 ) 
          # focal taxon去掉，不随机
          which.replace.nonfocal <- which.replace[!(which.replace %in% which.taxon)]
          #对样本内，不是0的，且不是focal taxon的物种进行随机化。
          perm.rel.d[j, which.replace.nonfocal] <- sample(rel.d[ j, which.replace.nonfocal]) 
        }
        
        # 算相关性
        cor.mat.null <- cor(perm.rel.d)
        
        # 保存每次对于focal物种迭代的结果
        perm.cor.vec.mat <- cbind(perm.cor.vec.mat, cor.mat.null[, which.taxon])
        
      }
      med.tax.cors <- cbind(med.tax.cors, apply(perm.cor.vec.mat, 1, median))
      # 运行情况
      if(which.taxon %% 20 == 0){print(which.taxon)}
    }
  }
}

#  实际相关-零模型相关
ifelse(use.custom.cors == T, {
  obs.exp.cors.mat <- custom.cor.mat.sub}, {
    obs.exp.cors.mat <- cor.mat.true - med.tax.cors
  }
)
diag(obs.exp.cors.mat) <- 0

#计算正负相关的连通性
connectedness.pos <- apply(obs.exp.cors.mat, 2, pos.mean)
connectedness.neg <- apply(obs.exp.cors.mat, 2, neg.mean)

# 根据定义，计算cohesion
cohesion.pos <- rel.d %*% connectedness.pos
cohesion.neg <- rel.d %*% connectedness.neg

####输出
output <- list(connectedness.neg, connectedness.pos, cohesion.neg, cohesion.pos)
names(output) <- c("Negative Connectedness", "Positive Connectedness", "Negative Cohesion", "Positive Cohesion")

print(output)

write.csv(connectedness.neg, file = "Negative Connectedness.csv")      
write.csv(connectedness.pos, file = "Positive Connectedness.csv")      
write.csv(cohesion.neg, file = "Negative Cohesion.csv")      
write.csv(cohesion.pos, file = "Positive Cohesion.csv")      
