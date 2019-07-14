#!/usr/bin/env Rscript

# http://stackoverflow.com/questions/13724063/if-else-constructs-inside-and-outside-functions
# https://stat.ethz.ch/R-manual/R-devel/library/base/html/sets.html

LibsRequired = c("VennDiagram", "grid", "graphics", "gridExtra", "gplots", "RColorBrewer", "gdata")

for (Lib in LibsRequired) {

	if (require(Lib, character.only = TRUE)) {
		print( paste(Lib, " is loaded correctly", sep="") ) 
	}else {
		stop( paste("Could not load ", Lib, ". Install library and re-run program.", sep="") )
	}

}

# running system command to know number of entries from CSV file
command<-paste("wc -l ", commandArgs()[7] , sep="")
x<-unlist(strsplit(system(command, intern = TRUE), "[[:space:]]", fixed = FALSE, perl = TRUE, useBytes = FALSE))
x<-as.numeric(x[1])

# Reading WithOutConcAvgDist_from_PC1_PC2_PC3
LUX_PC12 <- data.matrix(read.table(commandArgs()[7], header= FALSE,fill=TRUE,dec=".",sep="\t",row.names=1,col.names=c(0:x)))
	# converting triangular matrix to square matrix #upperTriangle(y)<-lowerTriangle(y) 
	upperTriangle(LUX_PC12) = 0
	LUX_PC12 = LUX_PC12 + t(LUX_PC12)
	colnames(LUX_PC12) = rownames(LUX_PC12)
LUX_PC12_dist <- as.dist(LUX_PC12)

# Reading WithOutConcAvgDist_from_PC1_PC2_PC3
LUX_PC123 <- data.matrix(read.table(commandArgs()[8], header= FALSE,fill=TRUE,dec=".",sep="\t",row.names=1,col.names=c(0:x)))
	# converting triangular matrix to square matrix #upperTriangle(y)<-lowerTriangle(y) 
	upperTriangle(LUX_PC123) = 0
	LUX_PC123 = LUX_PC123 + t(LUX_PC123)
	colnames(LUX_PC123) = rownames(LUX_PC123)
LUX_PC123_dist <- as.dist(LUX_PC123)

# Reading WithOutConcAvgDist_from_SimilarityScore
LUX_SS <- data.matrix(read.table(commandArgs()[9], header= FALSE,fill=TRUE,dec=".",sep="\t",row.names=1,col.names=c(0:x)))
	# converting triangular matrix to square matrix #upperTriangle(y)<-lowerTriangle(y) 
	upperTriangle(LUX_SS) = 0
	LUX_SS = LUX_SS + t(LUX_SS)
	colnames(LUX_SS) = rownames(LUX_SS)
LUX_SS_dist <- as.dist(LUX_SS)

# Saving dendrograms to pdfrow.names
OutFileName = paste ("./", commandArgs()[6], "/Dendrograms.pdf", sep="")
pdf(OutFileName, paper="a4")

# Dendrogram of LUX PC12
plot(hclust(LUX_PC12_dist, method = "complete", members=NULL), cex=1, ylab="a.u", main="LUX PC12")
# Heatmap of LUX PC12
my_palette = colorRampPalette(c("yellow", "blue"))(n=256)
heatmap.2(LUX_PC12, dendrogram="column", distfun = as.dist, hclustfun = hclust, key=TRUE, keysize=1, trace="none", density.info=c("none"), col= my_palette, symm=TRUE, cexRow=1, cexCol=1, cellnote=round(LUX_PC12, digits=2), main="LUX PC12" )

# Dendrogram of LUX PC123
plot(hclust(LUX_PC123_dist, method = "complete", members=NULL), cex=1, ylab="a.u", main="LUX PC123")
# Heatmap of LUX PC12
my_palette = colorRampPalette(c("yellow", "blue"))(n=256)
heatmap.2(LUX_PC123, dendrogram="column", distfun = as.dist, hclustfun = hclust, key=TRUE, keysize=1, trace="none", density.info=c("none"), col= my_palette, symm=TRUE, cexRow=1, cexCol=1, cellnote=round(LUX_PC123, digits=2), main="LUX PC123" )

# Dendrogram of LUX SS
plot(hclust(LUX_SS_dist, method = "complete", members=NULL), cex=1, ylab="a.u", main="LUX PC123")
# Heatmap of LUX SS
my_palette = colorRampPalette(c("yellow", "blue"))(n=256)
heatmap.2(LUX_SS, dendrogram="column", distfun = as.dist, hclustfun = hclust, key=TRUE, keysize=1, trace="none", density.info=c("none"), col= my_palette, symm=TRUE, cexRow=1, cexCol=1, cellnote=round(LUX_SS, digits=2), main="LUX SS" )


########################################################
# From Pear.R
########################################################

# Reading combined data file
DataFile = paste("./", commandArgs()[6], "/", commandArgs()[6], ".tab", sep="")
DataTable = data.frame(read.table(DataFile, header=TRUE, fill=FALSE, dec=".", sep="\t", row.names=1, stringsAsFactors = FALSE, strip.white = T))
len = length(colnames(DataTable))

# Collecting set-distance(s) as a data frame
SetDistanceTable <- matrix(0, nrow = len, ncol = len)
colnames(SetDistanceTable) <- colnames(DataTable)
row.names(SetDistanceTable) <- colnames(DataTable)

# Collecting set-similarity(s) as a data frame
SetSimilarityTable <- matrix(0, nrow = len, ncol = len)
colnames(SetSimilarityTable) <- colnames(DataTable)
row.names(SetSimilarityTable) <- colnames(DataTable)

count=1
lenAdelB_sum = 0
# cat("No. of Unique Lipids in all pairs of lipidomes\n")

i = 1
while ( i < len ){

	SetDistanceTable[i, i] <- 0.00
	SetSimilarityTable[i, i] <- 1

#	cat("No. of Lipids in",  (colnames(DataTable)[i]), ( length(which( !is.na(DataTable[,i]))) ), "\n", sep = " ")

	j = i+1
	while ( j <= length(colnames(DataTable)) ){

#		cat("No. of Lipids in",  (colnames(DataTable)[j]), ( length(which( !is.na(DataTable[,j]))) ), "\n", sep = " ")

		X = which( !is.na(DataTable[,i]))
		Y = which( !is.na(DataTable[,j]))

		lenA = length(which( !is.na(DataTable[,i])))
		lenB = length(which( !is.na(DataTable[,j])))

		lenAuB = length(union(X, Y))
		lenAnB = length(intersect(X, Y))

		lenSetDiffAB = length(setdiff(X, Y))
		lenSetDiffBA = length(setdiff(Y, X))

		lenAdelB = (lenSetDiffAB + lenSetDiffBA)

		result1 = round ( ( (lenSetDiffAB + lenSetDiffBA) / lenAuB ), digits = 2)
		result2 = round ( (lenAnB / lenAuB), digits = 2)
		result3 = round ( ( lenAnB / ( (lenA + lenB)/2 ) ), digits = 4  )
		result4 = (1-result3)

#		print(result3)
	
		SetDistanceTable[j, i] <- result4
		SetDistanceTable[i, j] <- result4

		SetSimilarityTable[j, i] <- result3
		SetSimilarityTable[i, j] <- result3

		j = j+1

#		cat(lenAdelB, " ", sep="")
		lenAdelB_sum = (lenAdelB_sum + lenAdelB)
		count=count+1
	}

	i=i+1
}

# print(SetSimilarityTable)

SetDistanceTable[i, i] <- 0.00
SetSimilarityTable[i, i] <- 1

# AvgAdelB = (lenAdelB_sum / count)
# cat ("\nSum \t No. of pairs \t Avg\n" )
# cat (lenAdelB_sum,  count, AvgAdelB, sep="\t" )
# cat("\n")

# Calculating Pearson's correlation coefficient between Samples based on Concentraions.
PearCor = cor (DataTable, method="pearson", use = "pairwise.complete.obs") # Pearson Correlation Coeffecient

# Calculating Pearson's correlation distance between Samples based on Concentraions.
PearCorDist = (function(x) x=(1-x)) (as.matrix(PearCor)) # Pearson Correlatin Distance = 1-r; r = Pearson Correlation Coeffecient

# Creating and plotting dendrogram based on Pearson Correlation Distance
text = bquote( italic(pearson~~distance) == 1-italic(pearson~~correlation)~(italic(r)) )
plot( hclust (as.dist(PearCorDist) , method="complete", members=NULL), main=text, ylab="a.u")

# Heatmap of Pearson Correlation distance
my_palette = colorRampPalette(c("yellow", "blue"))(n=256)
heatmap.2(PearCorDist, dendrogram="column", distfun = as.dist, hclustfun = hclust, key=TRUE, keysize=1, trace="none", density.info=c("none"), col= my_palette, symm=TRUE, cexRow=1, cexCol=1, cellnote=round(PearCorDist, digits=2), main="perarson distance" )

# Creating and plotting dendrogram based on set-similarity
text = bquote( italic(similarity~index) == over( italic(A)~intersect(italic(B)), {(italic(A)+italic(B))/2} ) )
plot( hclust (as.dist(SetSimilarityTable) , method="complete", members=NULL), main=text, ylab="a.u" )

# Heatmap of set-similairty NOT distance
heatmap.2(SetSimilarityTable, dendrogram="column", distfun = as.dist, hclustfun = hclust, key=TRUE, keysize=1, trace="none", density.info=c("none"), col=my_palette, symm=TRUE, cexRow=1, cexCol=1, cellnote=round(SetSimilarityTable, digits=2), main=text)

# Creating and plotting dendrogram based on set-distance
text = bquote( italic(distance~index) == 1~-~{over( italic(A)~intersect(italic(B)), {(italic(A)+italic(B))/2} ) } )
plot( hclust (as.dist(SetDistanceTable) , method="complete", members=NULL), main=text, ylab="a.u")

# Heatmap of set-distance scores
heatmap.2(SetDistanceTable, dendrogram="column", distfun = as.dist, hclustfun = hclust, key=TRUE, keysize=1, trace="none", density.info=c("none"), col=my_palette, symm=TRUE, cexRow=1, cexCol=1, cellnote=round(SetDistanceTable, digits=2), main=text)

dev.off()
