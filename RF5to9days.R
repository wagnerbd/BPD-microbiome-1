#read in full dataset
ra=read.table("C:/Users/brandie/Documents/Childrens Hospital/BPD/Gerber microbiome/Data/Rdata5to9days3.txt", header = TRUE, row.names = 1,sep = "\t")
matrix.dat  <- as.matrix(ra)
print(matrix.dat)
dim(ra)
print (ra[1:10, ])


##random forest
library(randomForest)
library(RColorBrewer)
set.seed(8085432)

help(randomForest)
##supervised
###################################################################
###BPD outcome - microbiome only - RA

seq.rf1 <- randomForest(BPD~ RA1 + RA2 + RA3 + RA4 + RA5 + RA6 + RA7 + RA8 + RA9 + RA10
 + RA11 + RA12 + RA13 + RA14 + RA15 + RA16 + RA17 + RA18  + TBL + ShannonH_Mean, ra, ntree=5000,  na.action=na.omit, proximity=TRUE)
print(seq.rf1)

varImpPlot(seq.rf1)
ralabels = c("Peptostreptococcus","Mycoplasma","Proteus", "Abiotrophia","Gemella","Klebsiella",
            "Fusobacterium","Neisseria","Enterococcus","Prevotella","Haemophilus","Corynebacterium",
            "Ureaplasma","Staphylococcus",
  "Escherichia-Shigella", "Rare" ,"Shannon diversity","TBL","Enterobacteriaceae","Streptococcus")
varImpPlot(seq.rf1,labels=ralabels, main="Variable Importance Plot")


summary(seq.rf1)
predict(seq.rf1)
importance(seq.rf1)

help(varImpPlot)

seq.r1f$votes
#plot(seq.rf1$votes[ ,2],group)

seq.rf1$type
FK <- cbind (genus, seq.rf1$importance)
seq.rf1$confusion
seq.rf1$forest$xbestsplit[,1]
seq.rf1$forest
 getTree(seq.rf, k=21, labelVar=TRUE)

MDSplot(seq.rf1, ra$BPD)
legend("topleft", legend=levels(ra$BPD), fill=brewer.pal(4, "Set1"), title="BPD status", bty="n") 


#plots the marginal effect of a variable on the class probability 
partialPlot(seq.rf1, ra, weight)
	## Looping over variables ranked by importance:
	imp <- importance(seq.rf1)
	impvar <- rownames(imp)[order(imp[, 1], decreasing=TRUE)]
	op <- par(mfrow=c(2, 3))
	for (i in c(1:6)) {
	partialPlot(seq.rf1, ra, impvar[i], xlab=impvar[i],
	main=" ")
	}
	par(op)

###################################################################
###BPD outcome - microbiome only - AA
set.seed(725)
seq.rf2 <- randomForest(BPD~ AA1 + AA2 + AA3 + AA4 + AA5 + AA6 + AA7 + AA8 + AA9 + AA10
 + AA11 + AA12 + AA13 + AA14 + AA15 + AA16 + AA17 + AA18  + ShannonH_Mean, ra, ntree=5000,  na.action=na.omit, proximity=TRUE)
print(seq.rf2)
varImpPlot(seq.rf2)

MDSplot(seq.rf2, ra$BPD)
legend("bottomright", legend=levels(ra$BPD), fill=brewer.pal(4, "Set1")) 

###################################################################

###BPD outcome - microbiome & clinical - RA

set.seed(256154)
seq.rf3 <- randomForest(BPD~ RA1 + RA2 + RA3 + RA4 + RA5 + RA6 + RA7 + RA8 + RA9 + RA10
 + RA11 + RA12 + RA13 + RA14 + RA15 + RA16 + RA17 + RA18 + TBL + ShannonH_Mean
 + gestational_age + zscore, ra, ntree=5000,  na.action=na.omit, proximity=TRUE)
print(seq.rf3)
varImpPlot(seq.rf3)

ralabels3 = c("Peptostreptococcus","Proteus","Mycoplasma", "Abiotrophia","Gemella","Klebsiella",
             "Fusobacterium","Neisseria","Enterococcus","Gestational age","Haemophilus","Prevotella","Corynebacterium",
             "Ureaplasma","Staphylococcus",
             "Escherichia-Shigella", "Zscore", "Rare" ,"Shannon diversity","TBL","Enterobacteriaceae","Streptococcus")
varImpPlot(seq.rf3,labels=ralabels3, main="Variable Importance Plot")



MDSplot(seq.rf3, ra$BPD)
legend("topleft", legend=levels(ra$BPD), fill=brewer.pal(4, "Set1"), title='BPD status', bty="n") 


seq.impute <- rfImpute(BPD ~ ., ra)
seq.rf3i <- randomForest(BPD~ RA1 + RA2 + RA3 + RA4 + RA5 + RA6 + RA7 + RA8 + RA9 + RA10
 + RA11 + RA12 + RA13 + RA14 + RA15 + RA16 + RA17 + RA18  + TBL + ShannonH_Mean
 + gestational_age + zscore + il_1b + il_2 + il_4 + il_5 + il_6 + il_8 + il_10 + il_12 + gm_csf
 + ifn_y + tnf_a + vegf, seq.impute, ntree=5000,  na.action=na.omit, proximity=TRUE)
print(seq.rf3i)
varImpPlot(seq.rf3i)
MDSplot(seq.rf3i, ra$BPD)
legend("bottomright", legend=levels(ra$BPD), fill=brewer.pal(4, "Set1")) 

###################################################################
###BPD outcome - microbiome & clinical - AA
set.seed(2014)
seq.rf4 <- randomForest(BPD~ AA1 + AA2 + AA3 + AA4 + AA5 + AA6 + AA7 + AA8 + AA9 + AA10
 + AA11 + AA12 + AA13 + AA14 + AA15 + AA16 + AA17 + AA18 + ShannonH_Mean
 + gestational_age + zscore + il_1b + il_2 + il_4 + il_5 + il_6 + il_8 + il_10 + il_12 + gm_csf
 + ifn_y + tnf_a + vegf, ra, ntree=5000,  na.action=na.omit, proximity=TRUE)
print(seq.rf4)
varImpPlot(seq.rf4)



###################################################################

###BPD outcome - microbiome & clinical REDUCED- RA

set.seed(320156)
seq.rf3r <- randomForest(BPD~ RA17 + RA3 + TBL + ShannonH_Mean + RA15 + zscore + RA5 +
        RA16 + RA18, ra, ntree=5000,  na.action=na.omit, proximity=TRUE)
print(seq.rf3r)
varImpPlot(seq.rf3r)




###################################################################
###PH outcome - microbiome & clinical - AA
set.seed(90565)
seq.rf5 <- randomForest(PH ~ AA1 + AA2 + AA3 + AA4 + AA5 + AA6 + AA7 + AA8 + AA9 + AA10
 + AA11 + AA12 + AA13 + AA14 + AA15 + AA16 + AA17 + AA18 + ShannonH_Mean
 + gestational_age + zscore + il_1b + il_2 + il_4 + il_5 + il_6 + il_8 + il_10 + il_12 + gm_csf
 + ifn_y + tnf_a + vegf, ra, ntree=5000,  na.action=na.omit, proximity=TRUE)
print(seq.rf5)
varImpPlot(seq.rf5)

MDSplot(seq.rf5, ra$PH)
legend("bottomright", legend=levels(ra$PH), fill=brewer.pal(4, "Set1")) 


###################################################################
###################################################################
##unsupervised
print(ra[1:10, 1:10])
dim(ra)
set.seed(743862)
seq2.rf <- randomForest( ~  AA1 + AA2 + AA3 + AA4 + AA5 + AA6 + AA7 + AA8 + AA9 + AA10
 + AA11 + AA12 + AA13 + AA14 + AA15 + AA16 + AA17 + AA18 + ShannonH_Mean
 + gestational_age + zscore + il_1b + il_2 + il_4 + il_5 + il_6 + il_8 + il_10 + il_12 + gm_csf
 + ifn_y + tnf_a + vegf, ra, ntree=10000,  na.action=na.omit)
print(seq2.rf)
varImpPlot(seq2.rf)

seq.mds <- MDSplot(seq2.rf, ra$BPD)
legend("topleft", legend=levels(ra$BPD), fill=brewer.pal(4, "Set1")) 

seq.mds <- MDSplot(seq2.rf, ra$state)
legend("topleft", legend=levels(ra$state), fill=brewer.pal(4, "Set1")) 

seq.mds <- MDSplot(seq2.rf, ra$PH)
legend("topleft", legend=levels(ra$PH), fill=brewer.pal(4, "Set1")) 


########################
dim(dat)
seq2.rf2 <- randomForest( ~  RA1 + RA2 + RA3 + RA4 + RA5 + RA6 + RA7 + RA8 + RA9 + RA10
 + RA11 + RA12 + RA13 + RA14 + RA15 + RA16 + RA17 + RA18 + TBL + ShannonH_Mean
 + gestational_age + zscore + il_1b + il_2 + il_4 + il_5 + il_6 + il_8 + il_10 + il_12 + gm_csf
 + ifn_y + tnf_a + vegf, ra, ntree=10000,  na.action=na.omit)
print(seq2.rf2)
varImpPlot(seq2.rf2)

seq.mds <- MDSplot(seq2.rf2, ra$BPD)
legend("topleft", legend=levels(ra$BPD), fill=brewer.pal(4, "Set1")) 

seq.mds <- MDSplot(seq2.rf2, ra$state)
legend("topleft", legend=levels(ra$state), fill=brewer.pal(4, "Set1")) 

seq.mds <- MDSplot(seq2.rf2, ra$PH)
legend("topleft", legend=levels(ra$PH), fill=brewer.pal(4, "Set1")) 

########################

seq2.rf3 <- randomForest( ~  RA1 + RA2 + RA3 + RA4 + RA5 + RA6 + RA7 + RA8 + RA9 + RA10
 + RA11 + RA12 + RA13 + RA14 + RA15 + RA16 + RA17 + RA18 + TBL + ShannonH_Mean, 
ra, ntree=10000,  na.action=na.omit, proximity=TRUE)
print(seq2.rf3)
varImpPlot(seq2.rf3)

seq.mds <- MDSplot(seq2.rf3, ra$BPD)
legend("bottomright", legend=levels(ra$BPD), fill=brewer.pal(4, "Set1")) 


seq.mds <- MDSplot(seq2.rf3, ra$state)
legend("topleft", legend=levels(ra$state), fill=brewer.pal(4, "Set1")) 

seq.mds <- MDSplot(seq2.rf3, ra$PH)
legend("topleft", legend=levels(ra$PH), fill=brewer.pal(4, "Set1")) 

library(cluster)
# calculate partitioning around the mediod from the unsupervised random forest
pamData <- pam(1 - seq2.rf3$proximity, k =3, diss = TRUE)
MDSplot(seq2.rf3, pamData$clustering, 
        bg=rainbow(length(levels(as.factor(pamData$clustering))))[unclass=pamData$clustering], 
        pch=c(21:25)[unclass=pamData$clustering], 
        palette=rep("black",5), 
        xlab="Dimension 1", 
        ylab="Dimension 2", 
        cex=.8)

# add a legend to the image
legend(
  "topleft",
  c("Cluster 1", "Cluster 2"), 
  pch=c(21:25), 
  col=rep("black", 3), 
  pt.bg=rainbow(3), 
  cex=.8)
########################
#microbiome only - non clr trans
seq.rf <- randomForest(grp~ ., sim2, ntree=10000,  na.action=na.omit, proximity=TRUE)



######################################
#unsupervised cytokine only
set.seed(32064321)

seq2.rf4 <- randomForest( ~  il_1b + il_2 + il_4 + il_5 + il_6 + il_8 + il_10 + il_12 + gm_csf
 + ifn_y + tnf_a + vegf, 
ra, ntree=10000,  na.action=na.roughfix, proximity=TRUE)
print(seq2.rf4 )
varImpPlot(seq2.rf4 )

seq.mds <- MDSplot(seq2.rf4 , ra$BPD)
legend("bottomright", legend=levels(ra$BPD), fill=brewer.pal(4, "Set1")) 


seq.mds <- MDSplot(seq2.rf4 , ra$state)
legend("topleft", legend=levels(ra$state), fill=brewer.pal(4, "Set1")) 

seq.mds <- MDSplot(seq2.rf4 , ra$PH)
legend("topleft", legend=levels(ra$PH), fill=brewer.pal(4, "Set1")) 

library(cluster)

# calculate partitioning around the mediod from the unsupervised random forest
pamData <- pam(1 - seq2.rf4$proximity, k =3, diss = TRUE)
MDSplot(seq2.rf4 , pamData$clustering, 
        bg=rainbow(length(levels(as.factor(pamData$clustering))))[unclass=pamData$clustering], 
        pch=c(21:25)[unclass=pamData$clustering], 
        palette=rep("black",5), 
        xlab="Dimension 1", 
        ylab="Dimension 2", 
        cex=.8)

# add a legend to the image
legend(
  "topright",
  c("Cluster 1", "Cluster 2", "Cluster 3"), 
  pch=c(21:25), 
  col=rep("black", 3), 
  pt.bg=rainbow(3), 
  cex=.8)


