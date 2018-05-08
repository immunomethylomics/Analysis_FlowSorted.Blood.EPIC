library(FlowSorted.Blood.EPIC)
library(minfi)
library(FlowSorted.Blood.450k)
library(EpiDISH)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

#For the original import using minfi and future import from GEO see FlowSorted.Blood.EPIC/scripts

library(ExperimentHub)
hub <- ExperimentHub()
query(hub, "FlowSorted.Blood.EPIC")
FlowSorted.Blood.EPIC <- hub[["EH1136"]]
FlowSorted.Blood.EPIC

data (IDOLOptimizedCpGs)
head (IDOLOptimizedCpGs)

RGsetTargets <- FlowSorted.Blood.EPIC[,
                                      FlowSorted.Blood.EPIC$CellType == "MIX"]

sampleNames(RGsetTargets) <- paste(RGsetTargets$CellType,
                                   seq_len(dim(RGsetTargets)[2]), sep = "_")
RGsetTargets

celltypes<-c("Bcell", "CD4T","CD8T","Mono","Neu", "NK")


pheno2<-as.data.frame(pData(FlowSorted.Blood.EPIC)) 

#General QC and supplementary figures
#Genetic heatmap control SNPs
library(pheatmap)
ann_colors = list(CellType=c( CD4T="#D95F02", CD8T ="#7570B3", Bcell="#1B9E77",NK="#E6AB02",
                              Neu="#66A61E", Mono="#E7298A", MIX="gray"),
                  Sex = c(F="red", M="blue"),
                  Ethnicity=c("African-American"='black', "East-Asian"="firebrick", "Indo-European"="lightgray", Mixed="green"))

#levels(pheno2$CellType)<-c( "CD4T" , "CD8T" , "Bcell", "NK", "Neu" , "Mono" , "MIX"   )
annotation_row<-data.frame(CellType=phen2o$CellType, Sex=pheno2$Sex, Ethnicity= pheno2$Ethnicity_wide,
                           #factor(pheno2$CellType, levels = c( "CD4T" , "CD8T" , "Bcell", "NK", "Neu" , "Mono" , "MIX"   ))
                           row.names = rownames(pheno2))
pheatmap(getSnpBeta(FlowSorted.Blood.EPIC), annotation_col = annotation_row, show_colnames = F, annotation_colors = ann_colors)

#Table one
library(tableone)
factorvars2<-c("CellType","Sex", "bmi_clas", "Ethnicity_wide", "Ethnic_self", "smoker")
design5<-as.data.frame(pheno2)
design5[factorvars2] <- lapply(design5[factorvars2], factor)

contvars<-c("Age","weight_kg", "height_m", "bmi", "purity")
design5[contvars] <- lapply(design5[contvars], as.numeric)

tableOne <- CreateTableOne(vars = c(factorvars2, contvars), data = design5)

print(tableOne,  quote = TRUE)


#Principal component regression analysis
cov<-data.frame(CellType=factor(pheno2$CellType),
    sex=factor(pheno2$Sex),
    slide=factor(pheno2$Slide),
    #array=factor(pheno2$Array),
    Ethnicity=pheno2$Ethnicity_wide,
    Age=pheno2$Age,
    weight=pheno2$weight_kg,
    height=pheno2$height_m,
    BMI=pheno2$bmi,
    Smoker=factor(pheno2$smoker),
    purity=pheno2$purity)
betas2<-preprocessRaw(FlowSorted.Blood.EPIC)
betas2<-getBeta(betas2)
betas3<-betas2[complete.cases(betas2),]
#Eliminate 0 variance
which(apply(betas3, 1, var)==0)#2 probes 0 variance
betas3<-betas3[apply(betas3, 1, var) != 0,]
dim(betas3)#[1] 866089     49 2 probes deleted
pcrplot(betas3, cov, npc=20)

npc=50
npc <- min(ncol(betas3), npc)
svd <- prcomp(t(betas3), center = TRUE, scale = TRUE, retx = TRUE,
              na.action = "na.omit")

screeplot(svd, npc, type = "barplot")

eigenvalue <- svd[["sdev"]]^2
prop <- (sum(eigenvalue[1:npc])/sum(eigenvalue)) * 100
cat("Top ", npc, " principal components can explain ", prop,
    "% of data \n    variation", "\n")
p <- ENmix:::lmmatrix(svd$x[, 1:npc], cov)
yaxis <- colnames(p)
plotp2<-function ( p, yaxis, xmax, title)
{
    par(mar = c(5, 4.2, 4, 2) + 0.1)
    plot(1, xlim = c(0, xmax), ylim = c(0, length(yaxis) + 1),
         type = "n", bty = "n", axes = FALSE, xlab = "Principal Component",
         ylab = "", main = title)
    axis(1, at = c(1:xmax), pos = 0.5, las = 1, lwd = 3)
    for (i in 1:length(yaxis)) {
        text(0.3, i, yaxis[i], xpd = TRUE, adj = 1)
    }
    for (i in 1:ncol(p)) {
        for (j in 1:nrow(p)) {
            pp <- p[j, i]
            colcode <- "white"
            if (pp <= 1e-09) {
                colcode = "darkred"
            }
            else if (pp <= 1e-04) {
                colcode = "red"
            }
            else if (pp <= 0.01) {
                colcode = "orange"
            }
            else if (pp <= 0.05) {
                colcode = "pink"
            }
            polygon(c(j - 0.5, j - 0.5, j + 0.5, j + 0.5), c(i -
                                                                 0.5, i + 0.5, i + 0.5, i - 0.5), col = colcode,
                    border = NA)
        }
    }
    legend("topright", c("<0.05", "<0.01", "<10E-5", "<10E-10"),
           col = c("pink", "orange", "red", "darkred"), pch = 15,
           pt.cex = 2, bty = "o", horiz = TRUE, xpd = TRUE)
}
plotp2(p, yaxis, 20, title = "Principal Component Regression Analysis")


#Table 1

anno450K<-getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
annoEPIC<-getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
table(IDOLOptimizedCpGs%in%rownames(anno450K))
table(IDOLOptimizedCpGs%in%rownames(annoEPIC))


#Purity supplementary figure

#install.packages("ggsignif")
#install.packages("ggpubr")
library(ggsignif)
library(ggpubr)
pheno<-pheno2
boxdat<-as.data.frame(pheno[pheno$CellType!="MIX",c("purity", "CellType")])
boxdat<-boxdat[order(boxdat$CellType),]
#boxdat$purity<-as.numeric(boxdat$purity)
my_comparisons<- list(c("Bcell", "Mono"),
                      c("CD4T", "Mono"),
                      c("CD8T", "Mono"),
                      c("Neu", "Mono"),
                      c("NK", "Mono"))#,
#c("CD8T", "NK"),
#c("Neu", "Bcell"))
#c("Bcell",  "CD4T", "CD8T", "Mono", "Neu", "NK")
ggboxplot(boxdat, x = "CellType", y = "purity",
          color = "CellType", palette = brewer.pal(8, "Dark2"),
          xlab = "Cell type", ylab = "Cell sorting purity %") +
    stat_compare_means(comparisons = my_comparisons, label.y = c(105, 103, 101, 102, 104))+
    stat_compare_means(label.y = 107) + theme(legend.position="none")


#Main figure 1 and information from table 1
#Heatmap comparison

anno<-getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

anno0<-as.data.frame(anno)

#Save the coefs for the heatmap.
#The fastest way is to debug(estimateCellCounts) and debug(estimateCellCounts2) the function and save the "coefs" object for the three steps
#Here we provide the probes from estimatecellcounts as RData objects for reproducibility
#load("/data/coefs450KQN.RData"))#600 probes from the 450kauto
#load("/data/coefsEPICauto.RData"))#600 probes from the EPICauto

counts450K <- estimateCellCounts(RGsetTargets, compositeCellType = "Blood",
                                 processMethod = "auto", probeSelect = "auto",
                                 cellTypes = c("CD8T","CD4T", "NK","Bcell","Mono","Neu"),
                                 referencePlatform = "IlluminaHumanMethylation450k",
                                 returnAll = TRUE, meanPlot = T, verbose = TRUE)


library(pheatmap)
#Load the coefs450KQN.RData
anno2<-anno0[rownames(anno0)%in% rownames(coefs450KQN),]
anno2$DNase<-ifelse(anno2$DNase_Hypersensitivity_NAME!="", "Yes", "No")
anno2$Enhancer<-ifelse(anno2$Phantom5_Enhancers!="", "Yes", "No")
annotation_row<-data.frame(DHS=anno2$DNase, Enhancer=anno2$Enhancer, row.names = rownames(anno2))
annotation_col<-data.frame(CellType=colnames(coefs450KQN), row.names = colnames(coefs450KQN))
table(anno2$DNase)
table(anno2$Enhancer)
table(anno2$Methyl450_Loci)
table(anno2$Relation_to_Island)
anno4<-strsplit(anno2$UCSC_RefGene_Group, ";")
#character 0 should be transformed to NA or the lapply will fail
#anno5a<-lapply(anno4, function(x) if(identical(x, character(0))) NA_character_ else x)
anno5a<-lapply(anno4, function(x) if(identical(x, character(0))) "Intergenic" else x)
anno5<-lapply(anno5a, `[[`, 1)
anno6<-factor(unlist(anno5), levels= c("TSS1500", "TSS200", "1stExon", "5'UTR","Body", "3'UTR", "Intergenic"))
table(anno6)


ann_colors = list(CellType=c(Bcell="#1B9E77", CD4T="#D95F02", CD8T ="#7570B3", Mono="#E7298A",
                             Neu="#66A61E", NK="#E6AB02"),
                  Enhancer=c(No="black", Yes="lightgray"))



pheatmap(coefs450KQN,color = colorRampPalette(c("yellow", "black", "blue"), space = "Lab")(128),
         annotation_row= annotation_row, annotation_col = annotation_col, fontsize = 14,
         show_rownames = F, show_colnames = T, legend = T, cluster_rows = T, annotation_colors = ann_colors,
         main = "600 probes automatic selection 450K (Reinius)")




#load("D:/Dropbox (Christensen Lab)/Deconvolution libraries/80 Dartmouth Deconvolution Libraries/IDOL Optimized Library for Deconvoluting Normal 6/FlowSorted.Blood.EPIC/coefsEPICauto.RData")


countsEPICauto<-estimateCellCounts2(RGsetTargets, compositeCellType = "Blood", 
                                    processMethod = "auto",
                                    probeSelect = "auto", 
                                    cellTypes = c("CD8T", "CD4T", "NK", "Bcell", 
                                                  "Mono", "Neu"), 
                                    referencePlatform = 
                                        "IlluminaHumanMethylationEPIC",
                                    referenceset = NULL,
                                    IDOLOptimizedCpGs =NULL, 
                                    returnAll = TRUE)
coefs<-coefsEPICauto
anno2<-anno[rownames(anno)%in% rownames(coefs),]
anno2$DNase<-ifelse(anno2$DNase_Hypersensitivity_NAME!="", "Yes", "No")
anno2$Enhancer<-ifelse(anno2$Phantom5_Enhancers!="", "Yes", "No")
annotation_row<-data.frame(DHS=anno2$DNase, Enhancer=anno2$Enhancer, row.names = rownames(anno2))
annotation_col<-data.frame(CellType=colnames(coefs), row.names = colnames(coefs))
table(anno2$DNase)
table(anno2$Enhancer)
table(anno2$Methyl450_Loci)
table(anno2$Relation_to_Island)
anno4<-strsplit(anno2$UCSC_RefGene_Group, ";")
#character 0 should be transformed to NA or the lapply will fail
#anno5a<-lapply(anno4, function(x) if(identical(x, character(0))) NA_character_ else x)
anno5a<-lapply(anno4, function(x) if(identical(x, character(0))) "Intergenic" else x)
anno5<-lapply(anno5a, `[[`, 1)
anno6<-factor(unlist(anno5), levels= c("TSS1500", "TSS200", "1stExon", "5'UTR","Body", "3'UTR", "Intergenic"))
table(anno6)

pheatmap(coefs,color = colorRampPalette(c("yellow", "black", "blue"), space = "Lab")(128),
         annotation_row= annotation_row, annotation_col = annotation_col, fontsize = 14,
         show_rownames = F, show_colnames = T, legend = T, cluster_rows = T, annotation_colors = ann_colors,
         main = "600 probes automatic selection EPIC")


#For the new method run first the automatic process
countsEPICIDOL<-estimateCellCounts2(RGsetTargets, compositeCellType = "Blood", 
                                    processMethod = "preprocessNoob",
                                    probeSelect = "IDOL", 
                                    cellTypes = c("CD8T", "CD4T", "NK", "Bcell", 
                                                  "Mono", "Neu"), 
                                    referencePlatform = 
                                        "IlluminaHumanMethylationEPIC",
                                    referenceset = NULL,
                                    IDOLOptimizedCpGs =IDOLOptimizedCpGs, 
                                    returnAll = TRUE)

anno2<-anno[rownames(anno)%in% IDOLOptimizedCpGs,]
anno2$DNase<-ifelse(anno2$DNase_Hypersensitivity_NAME!="", "Yes", "No")
anno2$Enhancer<-ifelse(anno2$Phantom5_Enhancers!="", "Yes", "No")
annotation_row<-data.frame(DHS=anno2$DNase, Enhancer=anno2$Enhancer, row.names = rownames(anno2))
annotation_col<-data.frame(CellType=colnames(countsEPICIDOL$compTable), row.names = IDOLOptimizedCpGs)
table(anno2$DNase)
table(anno2$Enhancer)
table(anno2$Methyl450_Loci)
table(anno2$Relation_to_Island)
anno4<-strsplit(anno2$UCSC_RefGene_Group, ";")
#character 0 should be transformed to NA or the lapply will fail
#anno5a<-lapply(anno4, function(x) if(identical(x, character(0))) NA_character_ else x)
anno5a<-lapply(anno4, function(x) if(identical(x, character(0))) "Intergenic" else x)
anno5<-lapply(anno5a, `[[`, 1)
anno6<-factor(unlist(anno5), levels= c("TSS1500", "TSS200", "1stExon", "5'UTR","Body", "3'UTR", "Intergenic"))
table(anno6)
library(FDb.InfiniumMethylation.hg19)
hm450 <- get450k()
#load("D:/OneDrive/OneDrive - Dartmouth College/Crossreactive/EPIC.manifest.rda")
#granges data is available from http://zwdzwd.github.io/InfiniumAnnotation
neartss<- getNearestTSS(EPIC.manifest)
neargene<-getNearestGene(EPIC.manifest)
neartss2<-neartss[rownames(neartss)%in%rownames(anno2),]
anno2$function_group<-anno6



pheatmap(countsEPICIDOL$compTable[IDOLOptimizedCpGs,3:8],color = colorRampPalette(c("yellow", "black", "blue"), space = "Lab")(128),
         annotation_row= annotation_row, annotation_col = annotation_col, fontsize = 14,
         show_rownames = F, show_colnames = T, legend = T, cluster_rows = T, annotation_colors = ann_colors,
         main = "DMR library IDOL selection EPIC")

library(gplots)
venn( list("450K (Reinius)"=rownames(coefs450KQN),"EPIC auto"=rownames(coefs),"EPIC IDOL"=IDOLOptimizedCpGs))



truevalues<-pheno[pheno$CellType == "MIX",celltypes] 




general_anno<-cbind.data.frame(countsEPICIDOL$compTable[order(row.names(countsEPICIDOL$compTable)),], anno2[order(rownames(anno2)),], neartss2[order(rownames(neartss2)),])

write.csv(general_anno, file="IDOL CpGs.csv")#Supplementary figure

#Comparison vs true values
#pheno<-as.data.frame(pData(RGsetTargets)) 
#pheno$clas<-paste("Method", substr(pheno$Subject.ID, 1,1))


pheno2$clas<-ifelse(pheno2$CellType=="MIX"& substr(pheno2$Sample_Name,1,1)=="A", "Method A",
                    ifelse(pheno2$CellType=="MIX"& substr(pheno2$Sample_Name,1,1)=="B", "Method B", as.character(pheno2$CellType)))

table(pheno$clas)


#Short plot
cellcount_EPICIDOL<-countsEPICIDOL$counts[,colnames(truevalues)]*100
rownames(cellcount_EPICIDOL)<-rownames(truevalues)#Use the original Sentrix ID

celltypes<-c( "CD4T","CD8T","Bcell", "NK","Neu","Mono")
#Remember the order for the colors look for the order if the results do not match the figure in the paper.
#Figure 1

par(mfrow=c(1,2))
barplot(t(truevalues[rownames(truevalues)%in% rownames(pheno2[pheno2$clas=="Method A",]),]), col= colorcelltype[celltypes], names.arg = seq(1:6), xlab = "Method A")
barplot(t(truevalues[rownames(truevalues)%in% rownames(pheno2[pheno2$clas=="Method B",]),]), col= colorcelltype[celltypes], names.arg = seq(1:6), legend.text = celltypes, xlab = "Method B")




par(mfrow=c(2,3))
for (i in celltypes) {
    #error <- truevalues[,i] - counts[,i]*100
    
    if (i=="Neu"){
        plot(truevalues[rownames(truevalues)%in% rownames(pheno2[pheno2$clas=="Method A",]),i],
             cellcount_EPICIDOL[rownames(cellcount_EPICIDOL)%in% rownames(pheno2[pheno2$clas=="Method A",]),i], main=paste(i),
             xlab="",ylab="",#paste("True", i, "%"), ylab=paste("Estimated", i, "%"),
             xlim=c(0,80), ylim=c(0,80),
             pch = 21,  col="black", bg=colorcelltype[i], #pch = 20,#as.numeric(mix),
             cex=2.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
        par(new=TRUE)
        plot(truevalues[rownames(truevalues)%in% rownames(pheno2[pheno2$clas=="Method B",]),i],
             cellcount_EPICIDOL[rownames(cellcount_EPICIDOL)%in% rownames(pheno2[pheno2$clas=="Method B",]),i], main=paste(i),
             xlab="",ylab="",#paste("True", i, "%"), ylab=paste("Estimated", i, "%"),
             xlim=c(0,80), ylim=c(0,80),
             pch = 22,  col="black", bg=colorcelltype[i], #pch = 20,#as.numeric(mix),
             cex=2.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
        lm.D9 <- lm(cellcount_EPICIDOL[,i]~truevalues[,i] )
        rmse(lm.D9$residuals)
        summary(lm.D9)$r.squared
        abline(lm.D9, col=colorcelltype[i])
        abline(a = 0, b=1, col="black", lwd=1, lty=2)
        text(10, 65, substitute(paste(R^2,m), list(m=round(summary(lm.D9)$r.squared[1]*100,1))), cex = 1.5)
        text(10, 55, substitute(paste("RMSE",m),list(m=round(rmse(lm.D9$residuals)[1],2))), cex = 1.5)
        
    }
    else {
        plot(truevalues[rownames(truevalues)%in% rownames(pheno2[pheno2$clas=="Method A",]),i],
             cellcount_EPICIDOL[rownames(cellcount_EPICIDOL)%in% rownames(pheno2[pheno2$clas=="Method A",]),i], main=paste(i),
             xlab="",ylab="",#paste("True", i, "%"), ylab=paste("Estimated", i, "%"),
             xlim=c(0,30), ylim=c(0,30),
             pch = 21,  col="black", bg=colorcelltype[i], #pch = 20,#as.numeric(mix),
             cex=2.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
        par(new=TRUE)
        plot(truevalues[rownames(truevalues)%in% rownames(pheno2[pheno2$clas=="Method B",]),i],
             cellcount_EPICIDOL[rownames(cellcount_EPICIDOL)%in% rownames(pheno2[pheno2$clas=="Method B",]),i], main=paste(i),
             xlab="",ylab="",#paste("True", i, "%"), ylab=paste("Estimated", i, "%"),
             xlim=c(0,30), ylim=c(0,30),
             pch = 22,  col="black", bg=colorcelltype[i], #pch = 20,#as.numeric(mix),
             cex=2.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
        lm.D9 <- lm(cellcount_EPICIDOL[,i]~truevalues[,i] )
        rmse(lm.D9$residuals)
        summary(lm.D9)$r.squared
        abline(lm.D9, col=colorcelltype[i])
        abline(a = 0, b=1, col="black", lwd=1, lty=2)
        text(5, 27, substitute(paste(R^2,m), list(m=round(summary(lm.D9)$r.squared[1]*100,1))), cex = 1.5)
        text(5, 23, substitute(paste("RMSE",m),list(m=round(rmse(lm.D9$residuals)[1],2))), cex = 1.5)
        
    }
}


#All the comparisons (supplementary figure information)
counts1<-counts450K$counts*100
cellcount_EPICauto<-countsEPICauto$counts*100
mix<-rep(1,12)
par(mfrow=c(2,2))
for (i in celltypes) {
    #error <- truevalues[,i] - counts[,i]*100
    
    
    
    if (i=="Neu"){
        plot(truevalues[,i], counts1[,i], main="Auto selection, CP/QP, Reinius 450K",
             xlab=paste("True", i, "%"), ylab=paste("Estimated", i, "%"), xlim=c(0,70), ylim=c(0,70),
             col=colorcelltype[i], pch = as.numeric(mix) , cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
        lm.D9 <- lm(counts1[,i]~truevalues[,i])
        rmse(lm.D9$residuals)
        summary(lm.D9)$r.squared
        abline(lm.D9, col=colorcelltype[i])
        abline(a = 0, b=1, col="black", lwd=1, lty=2)
        text(40, 20, substitute(paste(R^2,m), list(m=round(summary(lm.D9)$r.squared[1]*100,1))), cex = 1.5)
        text(40, 10, substitute(paste("RMSE",m),list(m=round(rmse(lm.D9$residuals)[1],2))), cex = 1.5)
        
    }
    else if (i=="Bcell" | i=="CD8T"){
        plot(truevalues[,i], counts1[,i], main="Auto selection, CP/QP, Reinius 450K",
             xlab=paste("True", i, "%"), ylab=paste("Estimated", i, "%"), xlim=c(0,30), ylim=c(0,30),
             col=colorcelltype[i], pch = as.numeric(mix), cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
        lm.D9 <- lm(counts1[,i]~truevalues[,i])
        rmse(lm.D9$residuals)
        summary(lm.D9)$r.squared
        abline(lm.D9, col=colorcelltype[i])
        abline(a = 0, b=1, col="black", lwd=1, lty=2)
        text(20, 13, substitute(paste(R^2,m), list(m=round(summary(lm.D9)$r.squared[1]*100,1))), cex = 1.5)
        text(20, 8, substitute(paste("RMSE",m),list(m=round(rmse(lm.D9$residuals)[1],2))), cex = 1.5)
        
    }
    else if (i=="CD4T"){
        plot(truevalues[,i], counts1[,i], main="Auto selection, CP/QP, Reinius 450K",
             xlab=paste("True", i, "%"), ylab=paste("Estimated", i, "%"), xlim=c(0,30), ylim=c(0,30),
             col=colorcelltype[i], pch = as.numeric(mix), cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
        lm.D9 <- lm(counts1[,i]~truevalues[,i])
        rmse(lm.D9$residuals)
        summary(lm.D9)$r.squared
        abline(lm.D9, col=colorcelltype[i])
        abline(a = 0, b=1, col="black", lwd=1, lty=2)
        text(15, 13, substitute(paste(R^2,m), list(m=round(summary(lm.D9)$r.squared[1]*100,1))), cex = 1.5)
        text(15, 8, substitute(paste("RMSE",m),list(m=round(rmse(lm.D9$residuals)[1],2))), cex = 1.5)
        
    }
    
    else{
        plot(truevalues[,i], counts1[,i], main="Auto selection, CP/QP, Reinius 450K",
             xlab=paste("True", i, "%"), ylab=paste("Estimated", i, "%"), xlim=c(0,30), ylim=c(0,30),
             col=colorcelltype[i], pch = as.numeric(mix), cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
        lm.D9 <- lm(counts1[,i]~truevalues[,i])
        rmse(lm.D9$residuals)
        summary(lm.D9)$r.squared
        abline(lm.D9, col=colorcelltype[i])
        abline(a = 0, b=1, col="black", lwd=1, lty=2)
        text(18, 13, substitute(paste(R^2,m), list(m=round(summary(lm.D9)$r.squared[1]*100,1))), cex = 1.5)
        text(18, 8, substitute(paste("RMSE",m),list(m=round(rmse(lm.D9$residuals)[1],2))), cex = 1.5)
        
    }
    
    
    
    
    if (i=="Neu"){
        plot(truevalues[,i], cellcount_EPICauto[,i], main="Auto selection, CP/QP, EPIC",
             xlab=paste("True", i, "%"), ylab=paste("Estimated", i, "%"), xlim=c(0,70), ylim=c(0,70),
             col=colorcelltype[i], pch = as.numeric(mix), cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
        lm.D9 <- lm(cellcount_EPICauto[,i]~truevalues[,i])
        rmse(lm.D9$residuals)
        summary(lm.D9)$r.squared
        abline(lm.D9, col=colorcelltype[i])
        abline(a = 0, b=1, col="black", lwd=1, lty=2)
        text(40, 20, substitute(paste(R^2,m), list(m=round(summary(lm.D9)$r.squared[1]*100,1))), cex = 1.5)
        text(40, 10, substitute(paste("RMSE",m),list(m=round(rmse(lm.D9$residuals)[1],2))), cex = 1.5)
        
    }
    else if (i=="Bcell" | i=="CD8T"){
        plot(truevalues[,i], cellcount_EPICauto[,i], main="Auto selection, CP/QP, EPIC",
             xlab=paste("True", i, "%"), ylab=paste("Estimated", i, "%"), xlim=c(0,30), ylim=c(0,30),
             col=colorcelltype[i], pch = as.numeric(mix), cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
        lm.D9 <- lm(cellcount_EPICauto[,i]~truevalues[,i])
        rmse(lm.D9$residuals)
        summary(lm.D9)$r.squared
        abline(lm.D9, col=colorcelltype[i])
        abline(a = 0, b=1, col="black", lwd=1, lty=2)
        text(20, 13, substitute(paste(R^2,m), list(m=round(summary(lm.D9)$r.squared[1]*100,1))), cex = 1.5)
        text(20, 8, substitute(paste("RMSE",m),list(m=round(rmse(lm.D9$residuals)[1],2))), cex = 1.5)
        
    }
    else if (i=="CD4T"){
        plot(truevalues[,i], cellcount_EPICauto[,i], main="Auto selection, CP/QP, EPIC",
             xlab=paste("True", i, "%"), ylab=paste("Estimated", i, "%"), xlim=c(0,30), ylim=c(0,30),
             col=colorcelltype[i], pch = as.numeric(mix), cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
        lm.D9 <- lm(cellcount_EPICauto[,i]~truevalues[,i])
        rmse(lm.D9$residuals)
        summary(lm.D9)$r.squared
        abline(lm.D9, col=colorcelltype[i])
        abline(a = 0, b=1, col="black", lwd=1, lty=2)
        text(15, 13, substitute(paste(R^2,m), list(m=round(summary(lm.D9)$r.squared[1]*100,1))), cex = 1.5)
        text(15, 8, substitute(paste("RMSE",m),list(m=round(rmse(lm.D9$residuals)[1],2))), cex = 1.5)
        
    }
    
    else{
        plot(truevalues[,i], cellcount_EPICauto[,i], main="Auto selection, CP/QP, EPIC",
             xlab=paste("True", i, "%"), ylab=paste("Estimated", i, "%"), xlim=c(0,30), ylim=c(0,30),
             col=colorcelltype[i], pch = as.numeric(mix), cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
        lm.D9 <- lm(cellcount_EPICauto[,i]~truevalues[,i])
        rmse(lm.D9$residuals)
        summary(lm.D9)$r.squared
        abline(lm.D9, col=colorcelltype[i])
        abline(a = 0, b=1, col="black", lwd=1, lty=2)
        text(18, 13, substitute(paste(R^2,m), list(m=round(summary(lm.D9)$r.squared[1]*100,1))), cex = 1.5)
        text(18, 8, substitute(paste("RMSE",m),list(m=round(rmse(lm.D9$residuals)[1],2))), cex = 1.5)
        
    }
    
    
    if (i=="Neu"){
        plot(truevalues[,i], cellcount_EPICIDOL[,i], main="IDOL, CP/QP, EPIC",
             xlab=paste("True", i, "%"), ylab=paste("Estimated", i, "%"), xlim=c(0,70), ylim=c(0,70),
             col=colorcelltype[i], pch = as.numeric(mix), cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
        lm.D9 <- lm(cellcount_EPICIDOL[,i]~truevalues[,i] )
        rmse(lm.D9$residuals)
        summary(lm.D9)$r.squared
        abline(lm.D9, col=colorcelltype[i])
        abline(a = 0, b=1, col="black", lwd=1, lty=2)
        text(40, 20, substitute(paste(R^2,m), list(m=round(summary(lm.D9)$r.squared[1]*100,1))), cex = 1.5)
        text(40, 10, substitute(paste("RMSE",m),list(m=round(rmse(lm.D9$residuals)[1],2))), cex = 1.5)
        
    }
    else if (i=="Bcell" | i=="CD8T"){
        plot(truevalues[,i], cellcount_EPICIDOL[,i], main="IDOL, CP/QP, EPIC",
             xlab=paste("True", i, "%"), ylab=paste("Estimated", i, "%"), xlim=c(0,30), ylim=c(0,30),
             col=colorcelltype[i], pch = as.numeric(mix), cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
        lm.D9 <- lm(cellcount_EPICIDOL[,i]~truevalues[,i] )
        rmse(lm.D9$residuals)
        summary(lm.D9)$r.squared
        abline(lm.D9, col=colorcelltype[i])
        abline(a = 0, b=1, col="black", lwd=1, lty=2)
        text(20, 13, substitute(paste(R^2,m), list(m=round(summary(lm.D9)$r.squared[1]*100,1))), cex = 1.5)
        text(20, 8, substitute(paste("RMSE",m),list(m=round(rmse(lm.D9$residuals)[1],2))), cex = 1.5)
        
    }
    else if (i=="CD4T"){
        plot(truevalues[,i], cellcount_EPICIDOL[,i], main="IDOL, CP/QP, EPIC",
             xlab=paste("True", i, "%"), ylab=paste("Estimated", i, "%"), xlim=c(0,30), ylim=c(0,30),
             col=colorcelltype[i], pch = as.numeric(mix), cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
        lm.D9 <- lm(cellcount_EPICIDOL[,i]~truevalues[,i] )
        rmse(lm.D9$residuals)
        summary(lm.D9)$r.squared
        abline(lm.D9, col=colorcelltype[i])
        abline(a = 0, b=1, col="black", lwd=1, lty=2)
        text(15, 13, substitute(paste(R^2,m), list(m=round(summary(lm.D9)$r.squared[1]*100,1))), cex = 1.5)
        text(15, 8, substitute(paste("RMSE",m),list(m=round(rmse(lm.D9$residuals)[1],2))), cex = 1.5)
        
    }
    
    else{
        plot(truevalues[,i], cellcount_EPICIDOL[,i], main="IDOL, CP/QP, EPIC",
             xlab=paste("True", i, "%"), ylab=paste("Estimated", i, "%"), xlim=c(0,30), ylim=c(0,30),
             col=colorcelltype[i], pch = as.numeric(mix), cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
        lm.D9 <- lm(cellcount_EPICIDOL[,i]~truevalues[,i] )
        rmse(lm.D9$residuals)
        summary(lm.D9)$r.squared
        abline(lm.D9, col=colorcelltype[i])
        abline(a = 0, b=1, col="black", lwd=1, lty=2)
        text(18, 13, substitute(paste(R^2,m), list(m=round(summary(lm.D9)$r.squared[1]*100,1))), cex = 1.5)
        text(18, 8, substitute(paste("RMSE",m),list(m=round(rmse(lm.D9$residuals)[1],2))), cex = 1.5)
        
    }
    plot(0,type='n',axes=FALSE,ann=FALSE)
}

#Figure 2
library(reshape2)
counts1<-counts1[,colnames(truevalues)]
difcounts1<-counts1-truevalues
#counts2<-as.data.frame(counts2)
cellcount_EPICauto<-cellcount_EPICauto[,colnames(truevalues)]
difEPICauto<-cellcount_EPICauto-truevalues
cellcount_EPICIDOL<-cellcount_EPICIDOL[,colnames(truevalues)]
difEPICIDOL<-cellcount_EPICIDOL-truevalues

difcounts1<-as.data.frame(difcounts1)
difcounts1$Method<-"CP/CQ Reinius"
difEPICauto<-as.data.frame(difEPICauto)
difEPICauto$Method<-"CP/CQ EPIC auto"
difEPICIDOL<-as.data.frame(difEPICIDOL)
difEPICIDOL$Method<-"CP/CQ EPIC IDOL"


difall<-rbind.data.frame(difcounts1, difEPICauto, difEPICIDOL)
#Do not work well
difcounts1<-difcounts1[, c("CD4T" , "CD8T" , "Bcell" ,"NK" ,   "Neu"  , "Mono", "Method")]
difEPICauto<-difEPICauto[, c("CD4T" , "CD8T" , "Bcell" ,"NK" ,   "Neu"  , "Mono", "Method")]
difEPICIDOL<-difEPICIDOL[, c("CD4T" , "CD8T" , "Bcell" ,"NK" ,   "Neu"  , "Mono", "Method")]

#difrelerr<-rbind.data.frame(df1=c(difcounts1[,-7]/truevalues, difcounts1$Method), df2=c(difEPICauto[,-7]/truevalues, difEPICauto$Method), df3=c(difEPICIDOL[,-7]/truevalues, difEPICauto$Method))

summary(difcounts1[,-7]/truevalues)
summary(difEPICauto[,-7]/truevalues)
summary(difEPICIDOL[,-7]/truevalues)
colorcelltype2<-c( Mono="#E7298A",
                   Neu="#66A61E", NK="#E6AB02", Bcell="#1B9E77",CD8T ="#7570B3", CD4T="#D95F02")
par(mfrow=c(1,1))
difall<-difall[, c("Mono","Neu","NK","Bcell", "CD8T", "CD4T", "Method")]


difmethod<-data.frame(Reinius=unlist(difall[which(difall$Method=="CP/CQ Reinius"), -7]),
                      EPIC_auto=unlist(difall[which(difall$Method=="CP/CQ EPIC auto"), -7]),
                      EPIC_IDOL=unlist(difall[which(difall$Method=="CP/CQ EPIC IDOL"), -7]))

t.test(difmethod$EPIC_IDOL, difmethod$Reinius)
t.test(difmethod$EPIC_IDOL, difmethod$EPIC_auto)
for (i in celltypes){
    print(i)
    print(bartlett.test(difall[, i]~difall$Method))
}

df.melt<-melt(difmethod)
p<-bartlett.test(df.melt$value~df.melt$variable)
p$p.value

df.m2<-melt(difall)
df.m2$color<-ifelse(df.m2$variable=="Bcell","#1B9E77",
                    ifelse(df.m2$variable=="CD4T","#D95F02",
                           ifelse(df.m2$variable=="CD8T","#7570B3",
                                  ifelse(df.m2$variable=="Mono","#E7298A",
                                         ifelse(df.m2$variable=="Neu","#66A61E", "#E6AB02")))))


library(RColorBrewer)
library(ggplot2)
library(gridExtra)
library(plotrix)

library(grid)
library(ggplot2)
library(gridBase)
pushViewport(viewport(layout = grid.layout(4, 1)))
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
#The trick seems to be to call plot.new before you set par, otherwise it's liable to get confused and not correctly honour the settings. You also need to set new = TRUE so a new page isn't started when you call plot.

#Create figure window and layout
plot.new()
grid.newpage()
pushViewport(viewport(layout = grid.layout(20, 1)))

#par(mfrow=c(2,1))
par(mar=c(5.1, 4.3, 4.1, 2.1)) #bottom, left, top, and right.
#Draw first plot
pushViewport(viewport(layout.pos.col = 1, layout.pos.row = c(1:12)))
par(fig = gridFIG(), mar=c(3.1, 4.3, 2.1, 2.1), new = TRUE)




boxplot(df.m2$value~df.m2$variable+factor(df.m2$Method, levels=c("CP/CQ Reinius", "CP/CQ EPIC auto", "CP/CQ EPIC IDOL" )) ,
        horizontal=T, col=colorcelltype2,
        border=c("blue", "blue","blue","blue","blue","blue",
                 "red", "red","red","red","red","red",
                 "black", "black","black","black","black","black"),
        las=1, yaxt='n')
abline(v=c(-5,5), col="black", lwd=1, lty=2)
abline(v=0, col="black", lwd=1, lty=1)
#abline(h=c(6.5,12.5), col="black", lwd=1, lty=1)
popViewport()


pushViewport(viewport(layout.pos.col = 1, layout.pos.row = c(12:20)))
par(fig = gridFIG(), mar=c(5.1, 4.3, 2.1, 2.1), new = TRUE)

boxplot(difmethod, horizontal=T, xlab="Estimated (%) - True (%)",
        las=1, boxwex=0.25, border=c("blue", "red","black"))
abline(v=c(-5,5), col="black", lwd=1, lty=2)
abline(v=0, col="black", lwd=1, lty=1)

abline(h=c(6.5,12.5), col="black", lwd=1, lty=1)
popViewport()

par(mar=c(5.1, 4.3, 4.1, 2.1))

#Styling was performed outside R

#Enrichment supplementary files
library(missMethyl)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(limma)
IDOL_CpGs <- read.csv("IDOL CpGs.csv", header=T, row.names = 1)
#View(IDOL_CpGs)
celltypes
#[1] "CD4T"  "CD8T"  "Bcell" "NK"    "Neu"   "Mono"
library(BAGS)
#gmt files are avilable from GSEA website please notice the version of the database for reproducibility
MSigDB_C2_Curated <- ReadGMT("D:\\OneDrive\\OneDrive - Dartmouth College\\GSEA/c2.all.v6.1.entrez.gmt")#Curated
MSigDB_C7_Curated <- ReadGMT("D:\\OneDrive\\OneDrive - Dartmouth College\\GSEA/c7.all.v6.1.entrez.gmt")#Immune


for(i in celltypes){
    gst <- gometh(sig.cpg=rownames(IDOL_CpGs[IDOL_CpGs$Cell_deconv==i,]), all.cpg=rownames(anno), collection="GO", array.type = "EPIC")
    print(paste("GO", i))
    print(topGO(gst))
    
    gsa.MSigDB <- gsameth(sig.cpg=rownames(IDOL_CpGs[IDOL_CpGs$Cell_deconv==i,]), all.cpg=rownames(anno), collection=MSigDB_C2_Curated, array.type = "EPIC")
    gsa.MSigDB <- as.data.frame(gsa.MSigDB[order(gsa.MSigDB[,4]),])
    sigMSigDB <- gsa.MSigDB[gsa.MSigDB$FDR<0.05,]
    print(paste("GSEA 2", i))
    print(sigMSigDB)
    
    gsa.MSigDB <- gsameth(sig.cpg=rownames(IDOL_CpGs[IDOL_CpGs$Cell_deconv==i,]), all.cpg=rownames(anno), collection=MSigDB_C7_Curated, array.type = "EPIC")
    gsa.MSigDB <- as.data.frame(gsa.MSigDB[order(gsa.MSigDB[,4]),])
    sigMSigDB <- gsa.MSigDB[gsa.MSigDB$FDR<0.05,]
    print(paste("GSEA 7", i))
    print(sigMSigDB)
}

gst <- gometh(sig.cpg=rownames(IDOL_CpGs), all.cpg=rownames(anno), collection="GO", array.type = "EPIC")
print(paste("GO", i))
print(topGO(gst))
sigGO <- gst[gst$FDR<0.05 & gst$DE>9 & gst$N<2001,]
sigGO<-sigGO[order(sigGO$FDR),]
write.csv(sigGO, file="GO enrichment.csv")

gsa.MSigDB <- gsameth(sig.cpg=rownames(IDOL_CpGs), all.cpg=rownames(anno), collection=MSigDB_C7_Curated, array.type = "EPIC")
gsa.MSigDB <- as.data.frame(gsa.MSigDB[order(gsa.MSigDB[,4]),])
sigMSigDB <- gsa.MSigDB[gsa.MSigDB$FDR<0.05 & gsa.MSigDB$DE>9 & gsa.MSigDB$N<2001,]
print(paste("GSEA 7", i))
print(sigMSigDB)
write.csv(sigMSigDB, file="GSEA immune enrichment.csv")


#Figure 6

cellcpgs<-c("cg25939861","cg04162316", "cg14047092", "cg03860768", "cg02647842", "cg22451300")

betascpgs<-getBeta(preprocessRaw(RGsetReferences))
cellcpgs<-cellcpgs[order(cellcpgs)]
betascpgs<-betascpgs[order(rownames(betascpgs)),]
betascpgs<-betascpgs[rownames(betascpgs)%in%cellcpgs,]
head(betascpgs)
phenoref<-pheno2[rownames(pheno2)%in%colnames(betascpgs),]
identical(rownames(phenoref), colnames(betascpgs))

boxplot(t(betascpgs)~phenoref$CellType)

cellcpgs<-c("cg04162316","cg25939861", "cg03860768", "cg14047092", "cg22451300", "cg02647842")
titlecpgs<-c(cg04162316="RPTOR", cg25939861="CD8A", cg03860768="BLK", cg14047092="CLASP1", cg22451300="NFIA", cg02647842="SLFN5")
colorcelltype<-colorcelltype[celltypes]
phenoref$CellType<-as.character(phenoref$CellType)
phenoref$CellType<-factor(phenoref$CellType, levels = celltypes)
par(mfrow=c(2,3))
for(i in cellcpgs){
    boxplot(betascpgs[i,]~phenoref$CellType, ylim=c(0,1), xlab=paste(i), ylab=expression(beta-value), main=paste(titlecpgs[i]),
            cex=2.5, cex.lab=1.3, cex.axis=1.3, cex.main=1.5, cex.sub=1.5,
            boxfill= colorcelltype)
}





#Longitudinal dataset
#Datasets are available  in GEO please see the manuscript for information.
#load the dataset first using the example provided in make dataset.R on the package hre the object is called as longdata
#setwd("D:/Dropbox (Christensen Lab)/Deconvolution libraries/80 Dartmouth Deconvolution Libraries/IDOL Optimized Library for Deconvoluting Normal 6/FlowSorted.Blood.EPIC")
#setwd("D:/Longitudinal Constitutive Samples_10-2017")
#sheet<-read.metharray.sheet(getwd())
#longData <- read.metharray.exp(targets = sheet,extended = F)

pheatmap(getSnpBeta(longData))#, annotation_col = annotation_row, show_colnames = F, annotation_colors = ann_colors)

predlong<-estimateCellCounts2(longData, referencePlatform ="IlluminaHumanMethylationEPIC", probeSelect = "IDOL", processMethod = "preprocessNoob")
predlong<-estimateCellCounts2(longData, compositeCellType = "Blood", 
                              processMethod = "preprocessNoob",
                              probeSelect = "IDOL", 
                              cellTypes = c("CD8T", "CD4T", "NK", "Bcell", 
                                            "Mono", "Neu"), 
                              referencePlatform = 
                                  "IlluminaHumanMethylationEPIC",
                              IDOLOptimizedCpGs =IDOLOptimizedCpGs, 
                              returnAll = FALSE)
cellcountlong<-predlong$counts

predlong450<-estimateCellCounts(longData, cellTypes = celltypes, returnAll = T)
predlong450<-estimateCellCounts2(longData, compositeCellType = "Blood", 
                                 processMethod = "preprocessQuantile",
                                 probeSelect = "auto", 
                                 cellTypes = c("CD8T", "CD4T", "NK", "Bcell", 
                                               "Mono", "Neu"), 
                                 referencePlatform = 
                                     "IlluminaHumanMethylation450k",
                                 #IDOLOptimizedCpGs =IDOLOptimizedCpGs, 
                                 returnAll = FALSE)
#save(longData, cellcountlong, file="Longitudinal data.RData")

dim(predlong$counts)

timeseq<-c(1,
           42,
           92,
           129,
           154,
           202,
           233,
           254,
           315,
           351,
           378,
           132#This is the corrected date
)

#Choose the comparison EPIC or 450K and run the plots
predlong$counts100<-predlong$counts*100
predlong$counts450K100<-predlong450$counts*100
countpred<-data.frame(predlong$counts100)
countpred$NLR<-countpred$Neu/(countpred$CD8T+countpred$CD4T+countpred$NK+countpred$Bcell)
countpred$LMR<-(countpred$CD8T+countpred$CD4T+countpred$NK+countpred$Bcell)/countpred$Mono
countpred$CD4_CD8<-countpred$CD4T/countpred$CD8T
countpred$CD8NK_Mono<-countpred$CD8T+countpred$NK/countpred$Mono
countpred$CD8_Bcell<-countpred$CD8T/countpred$Bcell
ratiosNLR<-c("NLR", "LMR", "CD4_CD8", "CD8NK_Mono", "CD8_Bcell")
colorratios<-c(NLR="blue", LMR="red", CD4_CD8="purple", CD8NK_Mono="black", CD8_Bcell="brown")
countpred2<-data.frame(predlong$counts450K100)
countpred2$NLR<-countpred2$Neu/(countpred2$CD8T+countpred2$CD4T+countpred2$NK+countpred2$Bcell)
countpred2$LMR<-(countpred2$CD8T+countpred2$CD4T+countpred2$NK+countpred2$Bcell)/countpred2$Mono
countpred2$CD4_CD8<-countpred2$CD4T/countpred2$CD8T
countpred2$CD4_CD8<-ifelse(countpred2$CD4_CD8<0, 0, countpred2$CD4_CD8)

countpred2$CD8NK_Mono<-countpred2$CD8T+countpred2$NK/countpred2$Mono
countpred2$CD8_Bcell<-countpred2$CD8T/countpred2$Bcell
countpred2$CD8_Bcell<-ifelse(is.infinite(countpred2$CD8_Bcell), 0, countpred2$CD8_Bcell)
#countpred2$CD8_Bcell<-ifelse(is.na(countpred2$CD8_Bcell), 0, countpred2$CD8_Bcell)

ratiosNLR<-c("NLR", "LMR", "CD4_CD8", "CD8NK_Mono", "CD8_Bcell")
colorratios<-c(NLR="blue", LMR="red", CD4_CD8="purple", CD8NK_Mono="black", CD8_Bcell="brown")

#Overlapped plot Figure 4 of the paper

#Just the plot
#
par(mfrow=c(2,1))
#par(mfrow=c(1,2))
#timeseq<-seq(1:12)

plot(timeseq[c(1:4, 12, 5:10)], predlong$counts100[c(1:4, 12, 5:10),celltypes[1]],type = "o",col = colorcelltype[celltypes[1]], xlab = "Day of measurement", ylab = "Estimated cell proportion (%)",
     #main = "Longitudinal cell type proportion estimation using the EPIC IDOL DMR library",
     ylim=c(0,80), cex=1.5, bty="n")
for(i in celltypes[2:6]){
    lines(timeseq[c(1:4, 12, 5:10)],predlong$counts100[c(1:4, 12, 5:10),i], type = "o", col = colorcelltype[i], cex=1.5)
}

for(i in celltypes[1:6]){
    lines(timeseq[c(1:4, 12, 5:10)],predlong$counts450K100[c(1:4, 12, 5:10),i], type = "l", col = colorcelltype[i], cex=1.5, lty=2)
}
legend(300,83, # places a legend at the appropriate place
       celltypes, # puts text in the legend
       lty=c(1,1), # gives the legend appropriate symbols (lines)
       lwd=c(2.5,2.5),col=colorcelltype) # gives the legend lines the correct color and width

plot(timeseq[c(1:4, 12, 5:10)], countpred[c(1:4, 12, 5:10),ratiosNLR[1]],type = "o",col = colorratios[ratiosNLR[1]], xlab = "Day of measurement", ylab = "Estimated cell ratios",
     #main = "Longitudinal cell ratio changes EPIC DMR library",
     ylim=c(0,25), cex=1.5, bty="n")
for(i in ratiosNLR[2:5]){
    lines(timeseq[c(1:4, 12, 5:10)], countpred[c(1:4, 12, 5:10),i], type = "o", col = colorratios[i], cex=1.5)
}
for(i in ratiosNLR[1:5]){
    lines(timeseq[c(1:4, 12, 5:10)], countpred2[c(1:4, 12, 5:10),i], type = "l", col = colorratios[i], cex=1.5, lty=2)
}
legend(300,26, # places a legend at the appropriate place
       c("Neu/Lymphocyte", "Lymphocyte/Mono", "CD4T/CD8T", "CD8T+NK/Mono", "CD8_Bcell"), # puts text in the legend
       lty=c(1,1), # gives the legend appropriate symbols (lines)
       lwd=c(2.5,2.5),col=colorratios) # gives the legend lines the correct color and width


#bartlet test
difall2<-rbind.data.frame(predlong$counts100[c(1:4, 12, 5:10),], predlong$counts450K100[c(1:4, 12, 5:10),])

difall2$Method<-NA
difall2$Method[1:11]<-"CP/CQ EPIC IDOL"
difall2$Method[11:22]<-"CP/CQ Reinius"


for (i in celltypes){
    print(i)
    print(t.test(difall2[, i]~difall2$Method))
    print(bartlett.test(difall2[, i]~difall2$Method))
}


#ratios
difall2<-rbind.data.frame(countpred[c(1:4, 12, 5:10),ratiosNLR], countpred2[c(1:4, 12, 5:10),ratiosNLR])

difall2$Method<-NA
difall2$Method[1:11]<-"CP/CQ EPIC IDOL"
difall2$Method[11:22]<-"CP/CQ Reinius"


for (i in ratiosNLR){
    print(i)
    print(t.test(difall2[, i]~difall2$Method))
    print(bartlett.test(difall2[, i]~difall2$Method))
}


#FACS validation

#FACS information is also available in GEO consult the specific datasets for information. Repeat for the three datasets the longitudinal dataset is already loaded
#Figure 5 from the manuscript
#EPIC dataset the object is called epicdata
setwd("D:/Dropbox (Christensen Lab)/FACS EPIC reference")
library(minfi)
library(ENmix)
sheet<-read.metharray.sheet(getwd())
#For ENmix it should be a extended RGset
epicData <- read.metharray.exp(targets = sheet,extended = TRUE)
names <- pData(epicData)$Sample_Name
groups <- pData(epicData)$CellType
densityBeanPlot(epicData, sampNames=names, sampGroups=groups)
densityPlot(epicData, sampGroups=groups)
densityBeanPlot(epicData, sampNames=names, sampGroups=groups)
qcReport(epicData)
plotCtrl(epicData,IDorder=NULL)
qcscore2<-QCinfo(epicData, detPthre=0.000001, nbthre=3, CpGthre=0.05, samplethre=0.05,outlier=TRUE, distplot=T)

#450K Koestler et al the object is called Mix450KData
setwd("D:/Dropbox (Christensen Lab)/Deconvolution libraries/80 Dartmouth Deconvolution Libraries/Reconstruction Mixtures 450")
sheet4<-read.metharray.sheet("D:/Dropbox (Christensen Lab)/Deconvolution libraries/80 Dartmouth Deconvolution Libraries/Reconstruction Mixtures 450")
Mix450KData <- read.metharray.exp(targets = sheet4,extended = TRUE, force=T)
mdsPlot(Mix450KData, sampGroups = sheet4$CellType, pal = c25, legendPos = "bottomright", legendNCol= 8)
qcscore450K<-QCinfo(Mix450KData, detPthre=0.05, nbthre=3, CpGthre=0.05, samplethre=0.05,outlier=TRUE, distplot=T)

#EPIC dataset
cell.Mix <- which(epicData$CellType == "WBC")
RGsetTargets <- epicData[, cell.Mix]
sampleNames(RGsetTargets) <- paste(RGsetTargets$CellType,
                                   seq(along = cell.Mix), sep = "_")
RGsetTargets
# Step 3: use your favorite package for deconvolution.
# Deconvolute a target data set consisting of EPIC DNA methylation 
# data profiled in blood, using your prefered method estimateCellCounts 
# (minfi), or similar.
# You can also read in the IDOL optimized DMR library based on the EPIC 
# array.  This object is nothing more than a vector of length 450 consisting 
# of the names of the IDOL optimized CpGs.  These CpGs are used as the 
# backbone for deconvolution and were selected because their methylation 
# signature differs across the six normal leukocyte subtypes.
data(IDOLOptimizedCpGs)

countsEPIC<-estimateCellCounts2(RGsetTargets, compositeCellType = "Blood", 
                                processMethod = "preprocessNoob",
                                probeSelect = "IDOL", 
                                cellTypes = c("CD8T", "CD4T", "NK", "Bcell", 
                                              "Mono", "Neu"), 
                                referencePlatform = 
                                    "IlluminaHumanMethylationEPIC",
                                IDOLOptimizedCpGs =IDOLOptimizedCpGs, 
                                returnAll = FALSE)



celltypes<-c( "CD4T","CD8T","Bcell", "NK","Neu","Mono")

pheno<-as.data.frame(pData(RGsetTargets))
pheno<-pheno[,celltypes]
pheno<-pheno*100
colorcelltype2<-c(Mono="#E7298A",
                  Neu="#66A61E", NK="#E6AB02", Bcell="#1B9E77",CD8T ="#7570B3", CD4T="#D95F02")

colorcelltype<-rev(colorcelltype2)
estimates<-as.data.frame(countsEPIC$counts)
estimates$Lymph<-estimates$CD8T+estimates$CD4T+estimates$NK+estimates$Bcell 
estimates<-estimates*100

#All in one
par(mfrow=c(1,1))
plot(pheno[,i],
     estimates[,i], #main=paste(i),
     xlab="",ylab="",#paste("True", i, "%"), ylab=paste("Estimated", i, "%"),
     xlim=c(0,80), ylim=c(0,80),
     pch = 22,  col="black", bg=colorcelltype[i], #pch = 20,#as.numeric(mix),
     cex=2.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
for (i in celltypes) {
    #error <- truevalues[,i] - counts[,i]*100
    par(new=TRUE)
    plot(pheno[,i],
         estimates[,i], #main=paste(i),
         xlab="",ylab="",#paste("True", i, "%"), ylab=paste("Estimated", i, "%"),
         xlim=c(0,80), ylim=c(0,80),
         pch = 22,  col="black", bg=colorcelltype[i], #pch = 20,#as.numeric(mix),
         cex=2.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
    #
    
}

rmse <- function(error)
{
    sqrt(mean(error^2))
}
lm.D9 <- lm(unlist(as.list(estimates[, celltypes]))~unlist(as.list(pheno[,celltypes])) )

rmse(lm.D9$residuals)
summary(lm.D9)$r.squared
abline(lm.D9, col="blue")
abline(a = 0, b=1, col="black", lwd=1, lty=2)
#text(10, 65, substitute(paste(R^2,m), list(m=round(summary(lm.D9)$r.squared[1]*100,1))), cex = 1.5)
#text(10, 55, substitute(paste("RMSE",m),list(m=round(rmse(lm.D9$residuals)[1],2))), cex = 1.5)

testdf<-matrix(data=NA, nrow = 6, ncol = 2, dimnames = list(celltypes, c("R2", "RMSE")))

for (i in celltypes){
    lm.D9 <- lm(estimates[,i]~pheno[,i] )
    rmse(lm.D9$residuals)
    summary(lm.D9)$r.squared
    testdf[i,1]<-round(summary(lm.D9)$r.squared[1]*100,1)
    testdf[i,2]<-round(rmse(lm.D9$residuals)[1],2)
}


library(plotrix)
addtable2plot(60 ,40,testdf,bty="o",display.rownames=TRUE,hlines=F,
              vlines=F,title="")


#450K
#EPIC dataset
cell.Mix <- which(Mix450KData$CellType == "WBC")
RGsetTargets <- Mix450KData[, cell.Mix]
sampleNames(RGsetTargets) <- paste(RGsetTargets$CellType,
                                   seq(along = cell.Mix), sep = "_")
RGsetTargets
# Step 3: use your favorite package for deconvolution.
# Deconvolute a target data set consisting of EPIC DNA methylation 
# data profiled in blood, using your prefered method estimateCellCounts 
# (minfi), or similar.
# You can also read in the IDOL optimized DMR library based on the EPIC 
# array.  This object is nothing more than a vector of length 450 consisting 
# of the names of the IDOL optimized CpGs.  These CpGs are used as the 
# backbone for deconvolution and were selected because their methylation 
# signature differs across the six normal leukocyte subtypes.
data(IDOLOptimizedCpGs)

countsEPIC<-estimateCellCounts2(RGsetTargets, compositeCellType = "Blood", 
                                processMethod = "preprocessNoob",
                                probeSelect = "IDOL", 
                                cellTypes = c("CD8T", "CD4T", "NK", "Bcell", 
                                              "Mono", "Neu"), 
                                referencePlatform = 
                                    "IlluminaHumanMethylationEPIC",
                                IDOLOptimizedCpGs =IDOLOptimizedCpGs450klegacy, 
                                returnAll = FALSE)



celltypes<-c( "CD4T","CD8T","Bcell", "NK","Neu","Mono")

pheno<-as.data.frame(pData(RGsetTargets))
pheno<-pheno[,celltypes]
pheno<-pheno*100
colorcelltype2<-c(Mono="#E7298A",
                  Neu="#66A61E", NK="#E6AB02", Bcell="#1B9E77",CD8T ="#7570B3", CD4T="#D95F02")

colorcelltype<-rev(colorcelltype2)
estimates<-as.data.frame(countsEPIC$counts)
estimates$Lymph<-estimates$CD8T+estimates$CD4T+estimates$NK+estimates$Bcell 
estimates<-estimates*100

#All in one
par(mfrow=c(1,1))
plot(pheno[,i],
     estimates[,i], #main=paste(i),
     xlab="",ylab="",#paste("True", i, "%"), ylab=paste("Estimated", i, "%"),
     xlim=c(0,80), ylim=c(0,80),
     pch = 22,  col="black", bg=colorcelltype[i], #pch = 20,#as.numeric(mix),
     cex=2.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
for (i in celltypes) {
    #error <- truevalues[,i] - counts[,i]*100
    par(new=TRUE)
    plot(pheno[,i],
         estimates[,i], #main=paste(i),
         xlab="",ylab="",#paste("True", i, "%"), ylab=paste("Estimated", i, "%"),
         xlim=c(0,80), ylim=c(0,80),
         pch = 22,  col="black", bg=colorcelltype[i], #pch = 20,#as.numeric(mix),
         cex=2.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
    #
    
}

rmse <- function(error)
{
    sqrt(mean(error^2))
}
lm.D9 <- lm(unlist(as.list(estimates[, celltypes]))~unlist(as.list(pheno[,celltypes])) )

rmse(lm.D9$residuals)
summary(lm.D9)$r.squared
abline(lm.D9, col="blue")
abline(a = 0, b=1, col="black", lwd=1, lty=2)
#text(10, 65, substitute(paste(R^2,m), list(m=round(summary(lm.D9)$r.squared[1]*100,1))), cex = 1.5)
#text(10, 55, substitute(paste("RMSE",m),list(m=round(rmse(lm.D9$residuals)[1],2))), cex = 1.5)

testdf<-matrix(data=NA, nrow = 6, ncol = 2, dimnames = list(celltypes, c("R2", "RMSE")))

for (i in celltypes){
    lm.D9 <- lm(estimates[,i]~pheno[,i] )
    rmse(lm.D9$residuals)
    summary(lm.D9)$r.squared
    testdf[i,1]<-round(summary(lm.D9)$r.squared[1]*100,1)
    testdf[i,2]<-round(rmse(lm.D9$residuals)[1],2)
}


library(plotrix)
addtable2plot(60 ,40,testdf,bty="o",display.rownames=TRUE,hlines=F,
              vlines=F,title="")


#Repeat for the longitudinal samples summarize the NK+Bcell first

