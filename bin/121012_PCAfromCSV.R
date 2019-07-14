#!/usr/bin/env Rscript

# Performing PCA

	# loading necessary PCA library
	library(gdata) 

	# running system command to know number of molecules from CSV file
	command<-paste("wc -l ", commandArgs()[6] , sep="")
	x<-unlist(strsplit(system(command, intern = TRUE), "[[:space:]]", fixed = FALSE, perl = TRUE, useBytes = FALSE))
	x<-as.numeric(x[1])

	# Reading CSV file
	y<-data.matrix(read.table(commandArgs()[6], header= FALSE,fill=TRUE,dec=".",sep=",",row.names=1,col.names=c(0:x)))

	# converting triangular matrix to square matrix #upperTriangle(y)<-lowerTriangle(y) 
	upperTriangle(y)<-0
	y<-y+t(y)

	# doing PCA analysis
	pca2<-princomp(y) 

	# Plotting the contribution of each component to the Variance
	myVector1 <- (pca2$sdev^2 / sum(pca2$sdev^2))*100
	pdf(file=commandArgs()[7])
	barplot(myVector1[1:10], ylab="% Variance")
	box()

	# Writing first three components to a seperate file for future use
	write.table(round(pca2$scores[,1:3],digits=4), file=commandArgs()[8]) 

# Calculating min and max score for each component in PCA (All graphs are plotted in same scale, useful for visual comparison)
xmin<-min(pca2$scores[,1])
xmax<-max(pca2$scores[,1])
ymin<-min(pca2$scores[,2])
ymax<-max(pca2$scores[,2])

cat(xmin, "\n")
cat(xmax, "\n")
cat(ymin, "\n")
cat(ymax, "\n")

dev.off()

# warnings()
