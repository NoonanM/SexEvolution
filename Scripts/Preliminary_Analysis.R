#Load the necessary packages
library(ape)
library(corrplot)
library(phytools)
library(ggtree)

data <- read.csv("Data/SextbookV4.csv")

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

#Generate the phylogenetically independent contrasts
for(i in 7:ncol(data)){  #22, 61, and 81 struggle with the seq(by = X) interval
  
  #Store the indexed trait
  DATA <- data[,i]
  
  if(is.numeric(DATA)){
    #Round
    DATA <- round(DATA,2)
    
    #Min and max values
    MIN <- min(na.omit(DATA))
    MAX <- max(na.omit(DATA))
    
    #Sequence of numbers that spans the range of the data
    VALS <- seq(MIN,
                MAX,
                by = 0.01)
    
    #Colours that span the range of the data
    COL_VEC <- viridis::viridis(length(VALS))
    
    #Match colours with the data (NAs are grey)
    COLS <- vector("character", length(DATA))
    for(j in 1:length(DATA)){
      KEEPER <- which(DATA[j] == VALS)
      if(length(KEEPER)>0){COLS[j] <- COL_VEC[KEEPER]} else {COLS[j] <- "grey80"}
    }
    
    #Generate the plot
    path <- paste("Figures/", names(data)[i],".png", sep = "")
    png(filename=path,
        width = 6, height = 10, units = "in",
        res = 300)
    
    plot(phylogeny,
         tip.color=COLS,
         cex = 0.42)
    title(names(data)[i])
    
    legend_image <- as.raster(matrix(rev(COL_VEC), ncol=1))
    
    par(mar=c(4, 2, 10, 12),                                  ## set margins
        new=TRUE)
    plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '')
    lbsq <- seq.int(0, 1, l=5)
    axis(4, at=lbsq, pos=0.1, labels=F, col=0, col.ticks=1, tck=-.05)
    
    
    
    text(x=0.4, y = seq(0,1,l=5), labels = seq(MIN,MAX,l=5),
         cex = 0.75)
    rasterImage(legend_image, 0, 0, 0.1,1)
    
    dev.off()
  }
  
  
  if(is.character(DATA)){
    
    DATA <- as.factor(DATA)
    FAKE <- as.numeric(DATA) #Just for arranging the colours easily
    
    
    #Sequence of numbers that spans the range of the data
    VALS <- unique(DATA)
    
    #Colours that span the range of the data
    COL_VEC <- viridis::viridis(length(VALS))
    
    #Match colours with the data (NAs are grey)
    COLS <- vector("character", length(DATA))
    for(j in 1:length(DATA)){
      KEEPER <- which(DATA[j] == VALS)
      if(length(KEEPER)>0){COLS[j] <- COL_VEC[KEEPER]} else {COLS[j] <- "grey80"}
    }
    
    #Colour by trait value but some extra work to handle NA values
    ncolors <- length(unique(DATA))
    
    path <- paste("Figures/", names(data)[i],".png", sep = "")
    png(filename=path,
        width = 6, height = 10, units = "in",
        res = 300)
    
    plot(phylogeny,
         tip.color=COLS,
         cex = 0.42)
    title(names(data)[i])
    
    legend_image <- as.raster(matrix(rev(COL_VEC), ncol=1))
    
    par(mar=c(4, 2, 10, 12),                                  ## set margins
        new=TRUE)
    plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '')
    lbsq <- seq.int(0, 1, l=ncolors)
    axis(4, at=lbsq, pos=0.1, labels=F, col=0, col.ticks=1, tck=-.05)
    text(x=0.4, y = seq(0,1,l=ncolors), labels = unique(DATA),cex = 0.75)
    rasterImage(legend_image, 0, 0, 0.1,1)
    
    dev.off()
  }
  
}






#Keep only the numeric data
keepers <- unlist(lapply(data, is.numeric))  
data_nums <- data[,keepers]

png(filename="Figures/Correlations.png",
    width = 6.86, height = 6, units = "in",
    res = 300)
corrplot::corrplot(cor(data_nums,
                       use = "pairwise.complete.obs"),
                   tl.cex = 0.5,
                   tl.col = "black",
                   type = "upper",
                   diag = FALSE)
dev.off()


