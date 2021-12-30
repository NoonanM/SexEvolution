#Load the necessary packages
library(ape)
library(corrplot)
library(phytools)
library(ggtree)

data <- read.csv("Data/SextbookV3.csv")

species <- data$Species
species <- gsub(" ", "_",species)
data$Species <- species

#Import the phylogeny
TREES <- read.nexus("Data/Phylogenies/output.nex")

#Calculate the consensus tree
phylogeny <- consensus.edges(TREES)

'%ni%' <- Negate('%in%')
species[which(species %ni% phylogeny$tip.label)]

#Then drop any species we don't have data for and line everything up
phylogeny <- drop.tip(phylogeny, c(which(phylogeny$tip.label %ni% data$Species)))
data <- data[match(phylogeny$tip.label, data$Species),]
phylogeny <- drop.tip(phylogeny, c(which(phylogeny$tip.label %ni% data$Species)))
data <- data[match(phylogeny$tip.label, data$Species),]

#Generate the phylogentically independent contrasts
for(i in 7:ncol(data)){
  DATA <- data[,i]
  BAD <- which(is.na(DATA))
  GOOD <- which(!is.na(DATA))
  DATA <- DATA[GOOD]
  
  if(length(DATA)>0){
    if(is.numeric(data[,i])){
      #Colour by trait value but some extra work to handle NA values
      ncolors <- length(unique(DATA))
      COLS <- viridis::viridis(ncolors)
      my_df <- data.frame(trait=DATA, color=COLS[cut(DATA,ncolors)])
      my_df_2 <- data.frame(trait=data[,i], species = data$Species)
      df <- merge(my_df_2, my_df, by = "trait", all = T)
      df <- df[!duplicated(df),]
      df$color[is.na(df$color)] <- "black"
      
      #Generate the plot
      path <- paste("Figures/", names(data)[i],".png", sep = "")
      png(filename=path,
          width = 6, height = 10, units = "in",
          res = 600)
      
      plot(phylogeny,
           tip.color=df$color,
           cex = 0.42)
      title(names(data)[i])
      
      legend_image <- as.raster(matrix(rev(COLS), ncol=1))
      
      par(mar=c(4, 2, 10, 12),                                  ## set margins
          new=TRUE)
      plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '')
      lbsq <- seq.int(0, 1, l=5)
      axis(4, at=lbsq, pos=0.1, labels=F, col=0, col.ticks=1, tck=-.05)
      
      
      
      text(x=0.4, y = seq(0,1,l=5), labels = seq(min(DATA),
                                                 max(DATA),l=5),
           cex = 0.75)
      rasterImage(legend_image, 0, 0, 0.1,1)
      
      dev.off()
    }
    
    if(is.character(data[,i])){
      
      DATA <- as.factor(DATA)
      FAKE <- as.numeric(DATA) #Just for arranging the colours easily
      
      #Colour by trait value but some extra work to handle NA values
      ncolors <- length(unique(DATA))
      COLS <- viridis::viridis(ncolors)
      my_df <- data.frame(trait=DATA, color=COLS[cut(FAKE,ncolors)])
      my_df_2 <- data.frame(trait=data[,i], species = data$Species)
      df <- merge(my_df_2, my_df, by = "trait", all = T)
      df <- df[!duplicated(df),]
      df$color[is.na(df$color)] <- "black"
      
      path <- paste("Figures/", names(data)[i],".png", sep = "")
      png(filename=path,
          width = 6, height = 10, units = "in",
          res = 600)
      
      plot(phylogeny,
           tip.color=df$color,
           cex = 0.42)
      title(names(data)[i])
      
      legend_image <- as.raster(matrix(rev(COLS), ncol=1))
      
      par(mar=c(4, 2, 10, 12),                                  ## set margins
          new=TRUE)
      plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '')
      lbsq <- seq.int(0, 1, l=ncolors)
      axis(4, at=lbsq, pos=0.1, labels=F, col=0, col.ticks=1, tck=-.05)
      text(x=0.4, y = seq(0,1,l=ncolors), labels = levels(DATA),cex = 0.75)
      rasterImage(legend_image, 0, 0, 0.1,1)
      
      dev.off()
    }
    
  } #Closes the if statement for the length of data > 0
}






#Keep only the numeric data
keepers <- unlist(lapply(data, is.numeric))  
data_nums <- data[,keepers]

png(filename="Figures/Correlations.png",
    width = 6.86, height = 6, units = "in",
    res = 600)
corrplot::corrplot(cor(data_nums,
                       use = "pairwise.complete.obs"),
                   tl.cex = 0.5,
                   tl.col = "black")
dev.off()


