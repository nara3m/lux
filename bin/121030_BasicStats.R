#!/usr/bin/env Rscript

y<-read.table(commandArgs()[6], header= FALSE,fill=TRUE,dec=".",sep=",");

z<-summary(y)

cat(z[1], "\n")
cat(z[2], "\n")
cat(z[3], "\n")
cat(z[4], "\n")
cat(z[5], "\n")
cat(z[6], "\n")

# warnings()

