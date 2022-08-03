###################################################################
### Use R to make PCAs from binary matrix of TE insertions
###Created by Jennifer Korstian
###################################################################
# Theese steps were run interactively on the HPCC at TTU 
# operates on a binary matrix rows= TEinsertion, columns= species, 1= insertion present, 0= insertion absent

###################################################################
### Load Required modules for R & dependencies
###################################################################
module load gnu7 R udunits geos gdal proj

install.packages("adegenet")
BiocManager::install("ggtree")
library(ggtree)
library(adegenet)
library(ape)
###################################################################
### Create pca scatter plot for TEs
###################################################################
# read in the binary matrix file
mat <- read.table(file="/lustre/work/jkorstia/MELT/ves_phylip_21apr21.matrix", sep=",", header=FALSE, row.names=1)

# create a new genlight object from the binary matrix
x=new("genlight", mat)

# run pca analysis of the genlight object and store as variable
pca_one <- glPca(x)
# must enter the number of axes manually (I never bothered to figure out how to automate this)
Select the number of axes: 2

# calculate eigenvalues for plot
eig.perc <- 100*pca_one$eig/sum(pca_one$eig)
# check that everything worked
head(eig.perc)
# 41.360110 16.231084 13.999613  8.219161  5.645878  4.372774 # sample eigen values

# create scatterplot based on PCA and write out to file
jpeg('/lustre/work/jkorstia/MELT/RplotsPCA_noeigen.jpg', width=5, height=5, units='in', res=300)
myCol <- colorplot(pca_one$scores,pca_one$scores, transp=TRUE, cex=4, xlab="PC 1 (XX%)", ylab="PC 2 (XX%)")
#text(pca_one$scores, labels=row.names(pca_one$scores))
#add.scatter.eig(pca_one$eig[1:11],2,1,2, posi="bottomleft", inset=0.1, ratio=0.3)
dev.off()

###################################################################
### Create neighbor joining tree for TEs
###################################################################

#infer nj tree from distance matrix
tre <- nj(dist(as.matrix(x)))

# reroot the tree & ladderize
tre2 <- root(tre, out=4)
tre2 <- ladderize(tre2)

# create plot of tree and write out to file
jpeg('/lustre/work/jkorstia/MELT/nj_tree_small_dots.jpg', width=6.5, height=4, units='in', res=400)
myCol <- colorplot(pca_one$scores,pca_one$scores, transp=TRUE)
par(mar = c(1.1, 2.1, 2.1, 1.1))
plot(tre2, type="unrooted", show.tip=TRUE, cex=.7)
tiplabels(pch=20, col=myCol, cex=1)
dev.off()
