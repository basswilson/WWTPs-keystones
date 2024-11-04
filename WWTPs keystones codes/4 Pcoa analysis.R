##################################################
## Project: WWTPs keystones
##  File name: Pcoa.R
## Date: 2024-11-04
## Author: Wilson
## R_Version: R version 4.4.1 (2024-06-14 ucrt)
##################################################


#Load packages
library(vegan) #Package required for distance calculations
library(ggplot2) #Plotting package
library(ggExtra)

#Read data, typically a data table where rows are sample names and columns are OTUs
otu_raw <- read.csv(file="all.csv",header=T,check.names=FALSE ,row.names=1)
#Transpose the data due to the requirements of the ordination functions
otu <- t(otu_raw)
otu <- as.data.frame(otu)

#Read the grouping file
group <- read.csv("metadata.csv",header=T)
group[group == "NA"] <- NA
colnames(group)[1] <- "samples"
# #Modify column names
# colnames(group) <- c("samples","group")

#Filter out NA rows based on different groups, using dissolved oxygen (DO) as an example
yenv<-group$DO
sampleid<-which(!is.na(yenv))  #Find the indices of non-missing values in yenv and store them in the variable sampleid
otu <- otu[sampleid,]

#Calculate Bray-Curtis distance
otu.distance <- vegdist(otu)
#Perform PCoA analysis
pcoa <- cmdscale(otu.distance,eig=TRUE)
pc12 <- pcoa$points[,1:2]
pc <- round(pcoa$eig/sum(pcoa$eig)*100,digits=2) #Percentage of explained variance

###Plotting###
#Convert pc12 from matrix to data.frame
pc12 <- as.data.frame(pc12)
#Add samples variable to pc12
pc12$samples <- row.names(pc12)
head(pc12)

#Plot
p <- ggplot(pc12,aes(x=V1, y=V2))+ #Specify data, X-axis, Y-axis
  geom_point(size=3)+ #Create point plot and set size
  theme_bw() #Theme
p

#Merge plotting data with grouping data
df <- merge(pc12,group,by="samples")

###############
#Perform adonis test
group_adonis <- df[,"DO"]
adonis_result <- adonis2(otu.distance ~ group_adonis, data = otu)

#View adonis results
adonis_result

R2 <- adonis_result$R2[1]
p_value <- adonis_result$'Pr(>F)'[1]

result_df <- data.frame(R2 = R2, p_value = p_value)
print(result_df)
write.csv(result_df,"DO.adonis.csv")

#Add adonis results to the plot
dune_adonis <- paste0("Adonis R2= ",formatC(R2, format = "f", digits = 4), "; p-value: ", p_value)
dune_adonis

#Plot

#################### For each group, you need to add a color, and the following code is similar.
color=c("#FF5733", "#00A896", "#FFC300", "#900C3F", "#FF8C42", 
        "#3F51B5", "#FF5722", "#06D6A0", "#FF6384", "#2196F3", 
        "#FFCC33", "#C2185B", "#00BFA5", "#FF4081", "#311B92", 
        "#FF8F00", "#39CCCC") #Color variable





p1 <-ggplot(data=df,aes(x=V1,y=V2,
                        color=DO))+ #Specify data, X-axis, Y-axis, color
  theme_bw()+ #Theme setting
  geom_point(size=1.8)+ #Create point plot and set size
  theme(panel.grid = element_blank())+
  geom_vline(xintercept = 0,lty="dashed")+
  geom_hline(yintercept = 0,lty="dashed")+ #Dashed lines in the plot
  #geom_text(aes(label=samples, y=V2+0.03,x=V1+0.03,  vjust=0),size=3.5)+ #Add labels for data points
  # guides(color=guide_legend(title=NULL))+ #Remove legend title
  labs(x=paste0("PC1 ",pc[1],"%"),
       y=paste0("PC2 ",pc[2],"%"),
       caption = dune_adonis,)+ #Change x, y axis titles to contribution percentages
  stat_ellipse(data=df,
               geom = "polygon",level=0.9,
               linetype = 2,size=0.5,
               aes(fill=DO),
               alpha=0.2,
               show.legend = T)+
  scale_color_manual(values = color) + #Set colors for points
  scale_fill_manual(values = c("#FF5733", "#00A896", "#FFC300", "#900C3F", "#FF8C42", "#3F51B5", "#FF5722", 
                               "#06D6A0", "#FF6384", "#2196F3", "#FFCC33", 
                               "#C2185B", "#00BFA5", "#FF4081", "#311B92", "#FF8F00", "#39CCCC"))+
  theme(axis.title.x=element_text(size=12), #Modify X-axis title text
        axis.title.y=element_text(size=12,angle=90), #Modify Y-axis title text
        axis.text.y=element_text(size=10), #Modify X-axis tick label text
        axis.text.x=element_text(size=10), #Modify Y-axis tick label text
        panel.grid=element_blank(), #Hide grid lines
        legend.position = "bottom") # Arrange legend horizontally
p1 <- ggMarginal(p1, type = "boxplot", margins = "both", groupColour = TRUE, groupFill = TRUE)

p1

ggsave("DO.png", width = 5, height = 4, plot = p1 )
ggsave("DO.svg", width = 5, height = 4, plot = p1 )

###############
###############
#Perform ANOVA test
model <- aov(V1 ~ DO, data = df)

#View ANOVA results
anova_result <- summary(model)

#Extract P-value
F <- anova_result[[1]]$`F`[1]
p_value <- anova_result[[1]]$`Pr(>F)`[1]
result_df <- data.frame(F = F, p_value = p_value)
#Output P-value
print(result_df)

write.csv(result_df, file = "DO.ANOVA.csv")


