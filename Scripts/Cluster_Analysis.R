library(cluster)
library(dendextend)
library(viridis)

data <- read.csv("Data/SextbookV4.csv")

#Remove duplicate data on WTD "White-tailed deer (South America)"
data <- data[-57,]

#Select only the numeric columns
data2 <- data[,unlist(lapply(data, is.numeric))]

row.names(data2) <- data$Species
dist <- daisy(data2,metric = "gower")

hc <- hclust(dist)

dend <- as.dendrogram(hc)

plot(dend,
     cex = 0.6)



png(filename="Figures/Cluster_Diet.png",
    width = 6.86, height = 5, units = "in",
    res = 600)

#Label colours
data3 <- data[match(labels(dend), data$Species),]
LABS <- as.numeric(as.factor(data3$Diet))
COLS <- pals::ocean.curl(length(unique(LABS)))[LABS]
labels_colors(dend) <- COLS
labels_cex(dend) <-0.3

plot(dend,
     cex = 0.2,
     axes=F)

legend("top",
       inset=.02,
       title="Diet",
       legend = c(levels(as.factor(data3$Diet))),
       fill=pals::ocean.curl(length(unique(LABS))),
       horiz=TRUE,
       cex=0.6,
       box.lty=0)

dev.off()



png(filename="Figures/Cluster_Sociality.png",
    width = 6.86, height = 5, units = "in",
    res = 600)

#Label colours
data3 <- data[match(labels(dend), data$Species),]
LABS <- as.numeric(as.factor(data3$solitary_group))
COLS <- pals::ocean.curl(length(unique(LABS)))[LABS]
labels_colors(dend) <- COLS
labels_cex(dend) <-0.3

plot(dend,
     cex = 0.2,
     axes=F)

legend("top",
       inset=.02,
       title="Sociality",
       legend = levels(as.factor(data$solitary_group)),
       fill=pals::ocean.curl(length(unique(LABS))),
       horiz=TRUE,
       cex=0.6,
       box.lty=0)

dev.off()



png(filename="Figures/Cluster_IUCN.png",
    width = 6.86, height = 5, units = "in",
    res = 600)

#Label colours
data3 <- data[match(labels(dend), data$Species),]
LABS <- as.numeric(as.factor(data3$Endangered_Status))
COLS <- pals::ocean.curl(length(unique(LABS)))[LABS]
labels_colors(dend) <- COLS
labels_cex(dend) <-0.3

plot(dend,
     cex = 0.2,
     axes=F)

legend("top",
       inset=.02,
       title="IUCN Status",
       legend = levels(as.factor(data$Endangered_Status)),
       fill=pals::ocean.curl(length(unique(LABS))),
       horiz=TRUE,
       cex=0.6,
       box.lty=0)

dev.off()




png(filename="Figures/Cluster_Family.png",
    width = 6.86, height = 5, units = "in",
    res = 600)

#Label colours
data3 <- data[match(labels(dend), data$Species),]
LABS <- as.numeric(as.factor(data3$Family))
COLS <- pals::ocean.curl(length(unique(LABS)))[LABS]
labels_colors(dend) <- COLS
labels_cex(dend) <-0.3

plot(dend,
     cex = 0.2,
     axes=F,
     main = "Coloured by Family")

dev.off()




#Run it again with only female traits
data <- read.csv("Data/female_traits.csv")

#Remove duplicate data on WTD "White-tailed deer (South America)"
data <- data[-36,]

#Select only the numeric columns
data2 <- data[,unlist(lapply(data, is.numeric))]

row.names(data2) <- data$Species
dist <- daisy(data2,metric = "gower")
hc <- hclust(dist)
dend <- as.dendrogram(hc)

plot(dend,
     cex = 0.6)



png(filename="Figures/Cluster_Female_Traits.png",
    width = 6.86, height = 5, units = "in",
    res = 600)
#Label colours
data3 <- data[match(labels(dend), data$Species),]
LABS <- as.numeric(as.factor(data3$Family))
COLS <- pals::ocean.curl(length(unique(LABS)))[LABS]
labels_colors(dend) <- COLS
labels_cex(dend) <-0.3
#Generate plot
plot(dend,
     cex = 0.2,
     axes=F,
     main = "Coloured by Family")
dev.off()




#Run it again with only female traits
data <- read.csv("Data/Male_traits.csv")

#Remove duplicate data on WTD "White-tailed deer (South America)"
data <- data[-36,]

#Select only the numeric columns
data2 <- data[,unlist(lapply(data, is.numeric))]

row.names(data2) <- data$Species
dist <- daisy(data2,metric = "gower")
dist[is.na(dist)] <- 0
hc <- hclust(dist)
dend <- as.dendrogram(hc)

plot(dend,
     cex = 0.6)



png(filename="Figures/Cluster_Male_Traits.png",
    width = 6.86, height = 5, units = "in",
    res = 600)
#Label colours
data3 <- data[match(labels(dend), data$Species),]
LABS <- as.numeric(as.factor(data3$Family))
COLS <- pals::ocean.curl(length(unique(LABS)))[LABS]
labels_colors(dend) <- COLS
labels_cex(dend) <-0.3
#Generate plot
plot(dend,
     cex = 0.2,
     axes=F,
     main = "Coloured by Family")
dev.off()

