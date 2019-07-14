#!/usr/bin/env Rscript

# cat(commandArgs()[6],"\n")
# cat(commandArgs()[7],"\n")
# cat(commandArgs()[8],"\n")
# cat(commandArgs()[9],"\n")
# cat(commandArgs()[10],"\n")
# cat(commandArgs()[11],"\n")
# cat(commandArgs()[12],"\n")
# cat(commandArgs()[13],"\n")
# cat(commandArgs()[14],"\n")
# cat(commandArgs()[6],"\n")

# loading SVG library
library(RSVGTipsDevice)

# Reading Concentrations File
ConcInput = read.table(commandArgs()[6], header=FALSE, fill=TRUE, sep="\t", dec=".", row.names=1)

# Modifying Molecule names in ConcInput file
rownames(ConcInput) = gsub(":|;| ", "_", rownames(ConcInput), ignore.case = FALSE, perl = TRUE, fixed = FALSE, useBytes = FALSE)

# Reading PCA File
PCAInput = read.table(commandArgs()[7], header=TRUE, fill=TRUE, sep=" ", dec=".", row.names=1)

# Collecting Coordinates for Molecules of ConcInput from PCAInput
# Custom Function to perform grep on a given molecule in PCAInput file
myFunc1 = function(MolName){
	search_string = paste ("\\b", MolName, "\\b", sep="")
	return(as.character(grep(search_string, rownames(PCAInput), ignore.case=FALSE, value=TRUE)));
}

# Preforming grep on each molecule of ConcInput using lapply
NewPCAInput = unlist(sapply(rownames(ConcInput), myFunc1),use.names=FALSE)
#cat(NewPCAInput,"\n")

# Reading Classes list file
ClassList = readLines(commandArgs()[8], n = -1L, ok = TRUE, warn = TRUE, encoding = "unknown")	# reading "classes_list" file; It contains short names for each class (PC/PA etc). These names are used in plot legend 

# Reading Centroid list file
CentroidList = read.table(commandArgs()[9], header=FALSE, fill=TRUE, sep="\t", dec=".", row.names=1) # This file contains common name for each lipid (a table of "lipid short name" and "centroid full name")

# Reading Xmin Xmax Ymin and Ymax values
xmin = as.numeric(commandArgs()[10])
xmax = as.numeric(commandArgs()[11])
ymin = as.numeric(commandArgs()[12])
ymax = as.numeric(commandArgs()[13])
HighestConc = as.numeric(commandArgs()[14])

# Picking Molecules that do not belong to any class and will be grouping them as 'Others'
myOthers = NewPCAInput
each_class = 1;
while(each_class <= length(ClassList)) {
	search_string = paste(ClassList[each_class],"_",sep="")
	myOthers = grep(search_string, myOthers, ignore.case=FALSE, fixed=TRUE, value=TRUE, invert=TRUE);
	each_class = each_class + 1;
}

# Plot 1 : WitOut Conc

# Naming output file
devSVGTips(file=paste(commandArgs()[6],"_1",".svg",sep=""), toolTipMode=2, onefile=FALSE, title="plot") # Try with TRUE (overlaid plots) or FALSE (separate plots)

# Defining a color set (for plot1)
MyColors = c("green", "red", "aquamarine", "blue", "blueviolet", "brown", "burlywood", "cadetblue", "chartreuse", "chocolate", "cornflowerblue", "darkgoldenrod", "darkgreen", "darkmagenta", "darkturquoise", "deeppink", "dimgrey", "forestgreen", "darkorange", "darkorchid", "darkred")

myColFunction = function (myCol) {
	myNewCol = rgb(col2rgb(myCol, alpha=FALSE)[1], col2rgb(myCol, alpha=FALSE)[2], col2rgb(myCol, alpha=FALSE)[3], 75, maxColorValue=255)
	return(myNewCol)
}

# Creating an empty scatter plot area... (Advanced Settings)
layout(matrix(c(1,1,1,1,2), 1, 5, byrow = TRUE))
par(mar=c(5,4,1,1)) # bottom, left, top, right margins
plot(NULL,NULL, xlab="Principal component 1",ylab="Principal component 2", xlim=c(xmin, xmax), ylim=c(ymin, ymax))
abline(v=0, h=0, col="gray", lty=1)

# Initializing values
each_class = 1; 
i = 1; 

myTransparentColSet= c()
myNewClassList = c()

while(each_class <= length(ClassList)) {

	myTransparentCol = myColFunction(MyColors[i])
	myTransparentColSet = c(myTransparentColSet, myTransparentCol)
	
	search_string = paste("^", ClassList[each_class],"_",sep="")

	k = grep(search_string,NewPCAInput,ignore.case=FALSE, value=TRUE);
#	cat(search_string, " ", length(k), "\n"); For Debugging like print 

	myNewClass = paste(ClassList[each_class], " [", length(k), "]", sep="") 
	myNewClassList = c(myNewClassList, myNewClass)

	sapply(k,function(k){
		setSVGShapeToolTip(desc1=k, desc2=CentroidList[k,1])
		points(PCAInput[k,1],PCAInput[k,2], pch=16, cex=2, col=myTransparentCol)})

	each_class = each_class+1;
	i=i+1;
}

# Plotting myOthers
	myTransparentCol = myColFunction(MyColors[i])
	myTransparentColSet = c(myTransparentColSet,myTransparentCol)

	sapply(myOthers,function(myOthers){
		setSVGShapeToolTip(desc1=myOthers, desc2=CentroidList[myOthers,1])
		points(PCAInput[myOthers,1],PCAInput[myOthers,2], pch=16, cex=2, col=myColFunction(MyColors[i]))})

myNewClass = paste("Others", " [", length(myOthers), "]", sep="") 
myNewClassList = c(myNewClassList, myNewClass)

o = c(rep(0.1, length(myNewClassList)))
p = rev(c(1:length(myNewClassList)))
par(mar=c(5,0,1,0))

plot(o, p, xlab="",ylab="",xlim=c(0,length(myNewClassList)), axes=FALSE, frame.plot=FALSE, pch=15, cex=3, col=myTransparentColSet,text(o,p,myNewClassList, pos=4, offset=1))

dev.off()
