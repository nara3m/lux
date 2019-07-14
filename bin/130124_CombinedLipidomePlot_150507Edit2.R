#!/usr/bin/env Rscript

# loading SVG library
library(RSVGTipsDevice)

# Reading PCA File
PCAInput = read.table(commandArgs()[7], header=TRUE, fill=TRUE, sep=" ", dec=".", row.names=1)

# Reading Classes list file
ClassList = readLines(commandArgs()[8], n = -1L, ok = TRUE, warn = TRUE, encoding = "unknown")	# reading "classes_list" file; It contains short names for each class (PC/PA etc). These names are used in plot legend 

# Finding Xmin Xmax Ymin and Ymax values
xmin = min(PCAInput[,1])
xmax = max(PCAInput[,1])
ymin = min(PCAInput[,2])
ymax = max(PCAInput[,2])
zmin = min(PCAInput[,3])
zmax = max(PCAInput[,3])

# Plotting code

# 2D plot

# Naming output file
	devSVGTips(file=commandArgs()[6], toolTipMode=2, onefile=FALSE, title="plot") # Try with TRUE (overlaid plots) or FALSE (separate plots)
	
	# Defining a color set (for plot1)
	MyColors = c("green", "red", "blue", "darkmagenta", "darkturquoise", "brown", "burlywood", "cadetblue", "blueviolet", "chartreuse", "chocolate", "cornflowerblue", "darkgoldenrod", "darkgreen", "deeppink", "forestgreen", "darkorange", "darkorchid", "darkred", "aquamarine")
	
	myColFunction = function (myCol) {
		myNewCol = rgb(col2rgb(myCol, alpha=FALSE)[1], col2rgb(myCol, alpha=FALSE)[2], col2rgb(myCol, alpha=FALSE)[3], 75, maxColorValue=255)
		return(myNewCol)
	}
	
	# Creating an empty scatter plot area... (Advanced Settings)
	# par(mar=c(5, 4, 4, 11)+0.1, omi=c(2,2,2,2)) # bottom, left, top, right margins
	
	layout(matrix(c(1,1,1,1,2), 1, 5, byrow = TRUE))
	par(mar=c(5,4,1,1))
	plot(NULL,NULL, xlab="Principal component 1",ylab="Principal component 2", xlim=c(xmin, xmax), ylim=c(ymin, ymax))
	abline(v=0, h=0, col="gray", lty=1)
	
	# Initializing values
	each_class = 1;
	i = 1; 
	
	myColSet= c()
	myNewClassList = c()
	
	LengthOfOtherMolecules=0
	
	while(each_class <= length(ClassList)) {
	
		SearchString = paste("^", ClassList[each_class],"_",sep="")
		SearchStringLength = nchar(SearchString)
		MoleculeIndex = grep(SearchString, substring(row.names(PCAInput), 1, SearchStringLength), ignore.case=FALSE, value=FALSE);
	
		if (length(MoleculeIndex) > 10){
			myTransparentCol = myColFunction(MyColors[i])
			myColSet = c(myColSet, MyColors[i])
			i=i+1;
	
			myNewClass = paste(ClassList[each_class], " [", length(MoleculeIndex), "]", sep="") 
			myNewClassList = c(myNewClassList, myNewClass)
	
		}
		else {
			myTransparentCol = myColFunction("grey")
			LengthOfOtherMolecules = LengthOfOtherMolecules + length(MoleculeIndex)
		}
	
		sapply(MoleculeIndex,function(MoleculeIndex){
			setSVGShapeToolTip(desc1=row.names(PCAInput)[MoleculeIndex])
			points(PCAInput[MoleculeIndex,1],PCAInput[MoleculeIndex,2], pch=16, cex=2, col=myTransparentCol)})
	
		each_class = each_class+1;
	}
	
	# Plotting MoleculesThatDontBelongToAnyClass
	
		# Finding Molecules that do not belong to any class and will be grouping them as 'Others'
		MoleculesThatDontBelongToAnyClass = seq(1, length(row.names(PCAInput)), by=1)
	
		each_class = 1;
		while(each_class <= length(ClassList)) {
			SearchString = paste("^", ClassList[each_class],"_",sep="")
			SearchStringLength = nchar(SearchString)
			MoleculesThatDontBelongToThisClass = grep(SearchString, substring(row.names(PCAInput), 1, SearchStringLength), ignore.case=FALSE, value=FALSE, invert=TRUE)
	
			MoleculesThatDontBelongToAnyClass = intersect(MoleculesThatDontBelongToAnyClass, MoleculesThatDontBelongToThisClass)
			each_class = each_class + 1;
		}
	
	sapply(MoleculesThatDontBelongToAnyClass,function(MoleculesThatDontBelongToAnyClass){
		setSVGShapeToolTip(desc1=row.names(PCAInput)[MoleculesThatDontBelongToAnyClass])
		points(PCAInput[MoleculesThatDontBelongToAnyClass,1],PCAInput[MoleculesThatDontBelongToAnyClass,2], pch=16, cex=2, col=myColFunction("grey"))})
	
	# par(xpd=TRUE)
	# legend(1.5,1,ClassList, fill=MyColors[1:11],cex=1, box.lwd=2)
	
	LengthOfOtherMolecules = LengthOfOtherMolecules + length(MoleculesThatDontBelongToAnyClass)
	myNewClass = paste("Others", " [", LengthOfOtherMolecules, "]", sep="") 
	myNewClassList = c(myNewClassList, myNewClass)
	
	o = c(rep(0.1, length(myNewClassList)))
	p = rev(c(1:length(myNewClassList)))
	myColSet = c(myColSet, "grey")
	par(mar=c(5,0,1,0))
	plot(o, p, xlab="",ylab="",xlim=c(0,length(myNewClassList)), axes=FALSE, frame.plot=FALSE, pch=15, cex=3, col=myColSet, text(o,p,myNewClassList, pos=4, offset=1))
	
	dev.off()

# 3D plot

library(scatterplot3d)

# Naming output file
	devSVGTips(file=commandArgs()[9], toolTipMode=2, onefile=FALSE, title="plot") # Try with TRUE (overlaid plots) or FALSE (separate plots)
	
	# Defining a color set (for plot1)
	MyColors = c("green", "red", "blue", "darkmagenta", "darkturquoise", "brown", "burlywood", "cadetblue", "blueviolet", "chartreuse", "chocolate", "cornflowerblue", "darkgoldenrod", "darkgreen", "deeppink", "forestgreen", "darkorange", "darkorchid", "darkred", "aquamarine")
	
	myColFunction = function (myCol) {
		myNewCol = rgb(col2rgb(myCol, alpha=FALSE)[1], col2rgb(myCol, alpha=FALSE)[2], col2rgb(myCol, alpha=FALSE)[3], 75, maxColorValue=255)
		return(myNewCol)
	}
	
	# Creating an empty scatter plot area... (Advanced Settings)
	# par(mar=c(5, 4, 4, 11)+0.1, omi=c(2,2,2,2)) # bottom, left, top, right margins
	
	layout(matrix(c(1,1,1,1,2), 1, 5, byrow = TRUE))
	par(mar=c(5,4,1,1))
	s3d<-scatterplot3d("PC 1", "PC 2","PC 3",xlim=c(xmin,xmax), ylim=c(ymin,ymax), zlim=c(zmin,zmax))
	
	# Initializing values
	each_class = 1;
	i = 1; 
	
	myColSet= c()
	myNewClassList = c()
	
	LengthOfOtherMolecules=0
	
	while(each_class <= length(ClassList)) {
	
		SearchString = paste("^", ClassList[each_class],"_",sep="")
		SearchStringLength = nchar(SearchString)
		MoleculeIndex = grep(SearchString, substring(row.names(PCAInput), 1, SearchStringLength), ignore.case=FALSE, value=FALSE);
	
		if (length(MoleculeIndex) > 10){
			myTransparentCol = myColFunction(MyColors[i])
			myColSet = c(myColSet, MyColors[i])
			i=i+1;
	
			myNewClass = paste(ClassList[each_class], " [", length(MoleculeIndex), "]", sep="") 
			myNewClassList = c(myNewClassList, myNewClass)
	
		}
		else {
			myTransparentCol = myColFunction("grey")
			LengthOfOtherMolecules = LengthOfOtherMolecules + length(MoleculeIndex)
		}
	
		sapply(MoleculeIndex,function(MoleculeIndex){
			setSVGShapeToolTip(desc1=row.names(PCAInput)[MoleculeIndex])
			s3d$points3d(PCAInput[MoleculeIndex,1], PCAInput[MoleculeIndex,2], PCAInput[MoleculeIndex,3], pch=16, cex=2, col=myTransparentCol)})
	
		each_class = each_class+1;
	}
	
	# Plotting MoleculesThatDontBelongToAnyClass
	
		# Finding Molecules that do not belong to any class and will be grouping them as 'Others'
		MoleculesThatDontBelongToAnyClass = seq(1, length(row.names(PCAInput)), by=1)
	
		each_class = 1;
		while(each_class <= length(ClassList)) {
			SearchString = paste("^", ClassList[each_class],"_",sep="")
			SearchStringLength = nchar(SearchString)
			MoleculesThatDontBelongToThisClass = grep(SearchString, substring(row.names(PCAInput), 1, SearchStringLength), ignore.case=FALSE, value=FALSE, invert=TRUE)
	
			MoleculesThatDontBelongToAnyClass = intersect(MoleculesThatDontBelongToAnyClass, MoleculesThatDontBelongToThisClass)
			each_class = each_class + 1;
		}
	
	sapply(MoleculesThatDontBelongToAnyClass,function(MoleculesThatDontBelongToAnyClass){
		setSVGShapeToolTip(desc1=row.names(PCAInput)[MoleculesThatDontBelongToAnyClass])
		s3d$points3d(PCAInput[MoleculesThatDontBelongToAnyClass,1], PCAInput[MoleculesThatDontBelongToAnyClass,2], PCAInput[MoleculesThatDontBelongToAnyClass,3], pch=16, cex=2, col=myColFunction("grey"))})
	
	# par(xpd=TRUE)
	# legend(1.5,1,ClassList, fill=MyColors[1:11],cex=1, box.lwd=2)
	
	LengthOfOtherMolecules = LengthOfOtherMolecules + length(MoleculesThatDontBelongToAnyClass)
	myNewClass = paste("Others", " [", LengthOfOtherMolecules, "]", sep="") 
	myNewClassList = c(myNewClassList, myNewClass)
	
	o = c(rep(0.1, length(myNewClassList)))
	p = rev(c(1:length(myNewClassList)))
	myColSet = c(myColSet, "grey")
	par(mar=c(5,0,1,0))
	plot(o, p, xlab="",ylab="",xlim=c(0,length(myNewClassList)), axes=FALSE, frame.plot=FALSE, pch=15, cex=3, col=myColSet, text(o,p,myNewClassList, pos=4, offset=1))
	
	# par(mar=c(5, 4, 4, 2)+0.1)
	
	dev.off()
	
	warnings()
