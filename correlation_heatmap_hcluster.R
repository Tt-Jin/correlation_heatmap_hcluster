library(tidyverse)
library(reshape2)
library(factoextra)
library(grid)
library(cowplot)

## read data
df <- read.csv("data/data.csv", row.names = 1, check.names = FALSE)
group <- read.csv("data/group.csv", check.names = FALSE)


## set color
mycol <- c("#A69900", "#8F1D00", "#056258", "#6A478F")


## calculate correlation and hclust
df_log <- log10(df + 0.000001)

corr <- cor(df_log, method = "pearson")
dist_matrix <- as.dist(1 - corr)
hclust_tree  <- hclust(dist_matrix, method = "complete")

## 得到上三角矩阵
get_lower_tri<-function(corr){
  corr[upper.tri(corr)] <- NA
  return(corr)
}
## 得到下三角矩阵
get_upper_tri <- function(corr){
  corr[lower.tri(corr)]<- NA
  return(corr)
}

## 重排序相关矩阵
corr <- corr[hclust_tree$order, hclust_tree$order]
upper_tri <- get_upper_tri(corr)


# hcluter tree
t <- fviz_dend(hclust_tree, 
               show_labels = FALSE, 
               lwd = 0.4, axes = FALSE,  
               ylab = "", main = "Cluster Dendrogram",
               ggtheme = theme_void()) +
  scale_x_discrete(expand = c(0,0)) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 14))
t


# heatmap
melted_corr <- melt(upper_tri, na.rm = TRUE)
p1 <- ggplot(data = melted_corr, aes(Var2, Var1))+
  geom_tile(aes(fill = value))+
  scale_fill_gradient2(low = "#BBDEEC", high = "#802520", 
                       breaks = c(seq(-0.2, 1.0, 0.2))) +
  scale_y_discrete(position = "right") +
  labs(x = NULL, y = NULL) +
  theme_minimal() + 
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(size = 10, color = "black", 
                                   angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 10, color = "black"),
        legend.text = element_text(size = 10))+
  coord_fixed() 
p1


# add group infomation
p2 <- p1 + 
  geom_point(data = group, 
             aes(x = `Bile acid`, y = -0.5, 
                 color = factor(order, levels = c("priConjSum", "primarySum", "secConjSum", "secondarySum"))),
             size = 2) +
  geom_point(data = group, 
             aes(y = `Bile acid`, x = 42.5, 
                 color = factor(order, levels = c("priConjSum", "primarySum", "secConjSum", "secondarySum"))),
             size = 2) + 
  scale_color_manual(values = mycol, na.translate=FALSE) +
  expand_limits(x = c(0, 43), y = c(-1, 42.5)) +
  guides(fill = guide_legend(title = "Correlation p", reverse = TRUE),
         color = guide_legend(title = NULL))

p2


# add rect
p3 <- p2 + geom_path(data=data.frame(x=c(10.5,16.5,16.5),
                                     y=c(10.5,10.5,16.5)),
                     aes(x=x,y=y), color="black", size=0.7)+
  geom_path(data=data.frame(x=c(22.5,32.5,32.5),
                            y=c(22.5,22.5,32.5)),
            aes(x=x,y=y), color="black", size=0.7) +
  geom_polygon(data=data.frame(x= c(0.5,0.5,41.5),
                               y = c(0.5,41.5,41.5)),
               aes(x=x,y=y), fill="white", color="white") +
  geom_polygon(data=data.frame(x=c(9.5,10.5,16.5,15.5),
                               y=c(11.5,10.5,16.5,17.5)),
               aes(x=x,y=y), fill="#E3E3E3") +
  geom_polygon(data=data.frame(x=c(21.5,22.5,32.5,31.5),
                               y=c(23.5,22.5,32.5,33.5)),
               aes(x=x,y=y), fill="#E3E3E3")

p3


# Rotate 45 degrees
p4 <- ggdraw()+
  draw_plot(ggplotify::as.ggplot(p3, angle = -45),
            width = 1, height = 1) 
p4 



# merge tree into heatmap
pdf("correlation_heatmap_hcluster.pdf", width = 14, height = 8)
p4
vie1 <- viewport(width=0.56, height=0.15, x=0.46, y=0.795)
print(t, vp=vie1)
dev.off()




