#### LOAD LIBRARIES ####
library(tidyverse)
library(ggfortify)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(edgeR)
library(mdatools)

#### put all your data files in one same document path
file_path <- '/Users/sharp/Desktop/Yasmeen/raw data/'

#### SELECT SPECIES AND FILTER OPTIONS ####
species<-"TBMU"   # set as "TBMU" or "BLGU"
looseFilter<-0
stringentFilter<-0.9
showBatch<-FALSE

#### SPECIES-SPECIFIC FILE NAMES AND INFO ####
switch(species,
       TBMU={
         chemFile<-"chemical-TBMU.csv"
         pcrFile<-"PCR-TBMU.csv"
         site<-FALSE
         batch<-showBatch
         batchData<-c(
           rep("batch 1", 17),
           rep("batch 2", 13)
         )
         missingChem<-c(28)
       },
       BLGU={
         chemFile<-"chemical-BLGU.csv"
         pcrFile<-"PCR-BLGU.csv"
         site<-TRUE
         batch<-showBatch
         siteData<-c(
           rep("Qaqulluit NWA", 16),
           rep("Akpait NWA", 14)
         )
         batchData<-c(
           rep("batch 1", 18),
           rep("batch 2", 12)
         )
         missingChem<-c(9,23)
       },
       stop("not a recognized option")
)

#### PLOT AND ANNOTATION COLORS ####
parentPalette<-brewer.pal(n = 4, name = "Oranges")
alkylPalette<-brewer.pal(n = 4, name = "Purples")

ann_colors = list(Site=c("Akpait NWA"='#ff9e91',"Qaqulluit NWA"="#6bfff8"),
                  Group=c(PLMW=parentPalette[2],
                          PHMW=parentPalette[3],
                          PHET=parentPalette[4],
                          ALMW=alkylPalette[2],
                          AHMW=alkylPalette[3],
                          AHET=alkylPalette[4],
                          'Trace Elements'="#6b5145"),
                  Category=c('Trace Elements'='black',
                             'Parent PACs'='grey50',
                             'Alkylated PACs'='grey80'),
                  Batch=c("batch 1"='green',"batch 2"="red"))

#### FUNCTIONS ####
# function mdl_filter(): only retains rows that have > indicated % of values above MDL
mdl_filter<-function(data, mdl, min_prop){   
  highDetect<-sapply(rownames(data),function(x){
    sum(data[x,]>mdl[x,])/ncol(data)>min_prop
  })
  filterData<-data[highDetect,]
  return(filterData)
}

#### LOAD CHEMICAL DATA ####
chemRaw <- read.csv(paste0(file_path,chemFile), row.names = 'Chemical')
group <- read.csv(paste0(file_path,'Chemical grouping.csv'))
mdl <- read.csv(paste0(file_path,'mdl.csv'),row.names='chemical')
plot(rowSums(chemRaw==0))

#replace anything < MDL values with MDL
correctedData<-chemRaw
for(i in rownames(correctedData)){
  correctedData[i,correctedData[i,]<mdl[i,]]<-mdl[i,]
}

#### CHEMICAL HEATMAP ####
# load chemical data and chemical groupings
heatData<-mdl_filter(correctedData, mdl, looseFilter)
heatData<-log10(heatData)
colnames(heatData)<-sapply(strsplit(colnames(heatData),'\\.'),function(x)x[2])

heatGroup<-group[group$Compound%in%rownames(heatData),]
rownames(heatGroup)<-heatGroup$Compound

row_order<-c(match(rownames(heatGroup), rownames(heatData)))
heatData<-heatData[row_order,]

# row annotations
Group<-factor(heatGroup$Grouping, levels=unique(heatGroup$Grouping))
Category<-factor(heatGroup$Category, levels=unique(heatGroup$Category))

ha_row<-rowAnnotation(Group=Group,
                      Category=Category,
                      col=ann_colors[2:3],
                      annotation_name_side="top")

# column annotations (if site=TRUE, or batch=TRUE)
if(site){
  heatSite<-siteData[-missingChem]
  names(heatSite)<-colnames(heatData)
  ha_col<-HeatmapAnnotation(Site=heatSite,
                            col=ann_colors[1],
                            annotation_name_side="left")
}else{
  ha_col<-NULL
}

if(batch){
  heatBatch<-batchData[-missingChem]
  names(heatBatch)<-colnames(heatData)
  ha_col<-HeatmapAnnotation(
    Batch=heatBatch,
    col=ann_colors[4],
    annotation_name_side="left")
}

if(site&batch){
  ha_col<-HeatmapAnnotation(
    Batch=heatBatch,
    Site=heatSite,
    col=ann_colors[c(1,4)],
    annotation_name_side="left")
}



# concentration colour palette
color_fun<- colorRamp2(c(min(heatData),0,max(heatData)), c("cornflowerblue", "#EEEEEE", "red"))

#plot heatmap
ht1<-Heatmap(as.matrix(heatData),name = "log10 (ng/g)",
             top_annotation = ha_col,
             right_annotation = ha_row,
             show_row_dend=F,
             clustering_distance_columns = "manhattan",
             column_title = paste0(species, ' Chemical Residues'),
             cluster_rows = F,
             column_names_gp = gpar(fontsize = 10),
             row_names_gp = gpar(fontsize = 8),
             col = color_fun)
ht1

#### CHEMICAL PCA ####
# filter and transcform chemical data
pcaData<-mdl_filter(correctedData, mdl, looseFilter)
pcaData<-t(log10(pcaData))
pcaData<-data.frame(pcaData, check.names =FALSE)
rownames(pcaData)<-sapply(strsplit(rownames(pcaData),'\\.'),function(x)x[2])

# pca model
pcaMod<-prcomp(pcaData,center=TRUE,scale=TRUE)

# top loadings
top5_pc1<-names(pcaMod$rotation[order(abs(pcaMod$rotation[,1]), decreasing=TRUE),1])[1:5]
top5_pc2<-names(pcaMod$rotation[order(abs(pcaMod$rotation[,2]), decreasing=TRUE),2])[1:5]
top5<-unique(c(top5_pc1, top5_pc2))

pcaMod2<-pcaMod
pcaMod2$rotation<-pcaMod2$rotation[top5,]
pcaMod2$center<-pcaMod2$center[top5]
pcaMod2$scale<-pcaMod2$scale[top5]


# pca color and frame options
pcaCol<-ifelse(site|batch, "legend", "red")
pcaFrame<-ifelse(site, TRUE, ifelse(batch, TRUE , FALSE))  

if(site){
  pcaData$legend<-siteData[-missingChem]
}

if(batch){
  pcaData$legend<-batchData[-missingChem]
}

# plot pca
p<-autoplot(pcaMod2, 
            # optional arguments for colony of origin info
            colour=pcaCol,
            frame=pcaFrame,
            # remaining common arguments
            shape = 16, size = 4,
            data=pcaData,label=F, 
            loadings = TRUE, loadings.label = TRUE,
            loadings.colour = 'darkgrey',
            loadings.label.colour = 'darkgrey',
            loadings.label.size = 4,
            main=paste0('PCA of ',species,' chemical burdens')) + 
  geom_text(aes(label = rownames(pcaData)), size = 4,
            nudge_y = rep(0.025,dim(pcaData)[1])) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black",fill=NA))
p


#### LOAD GENE EXPRESSION DATA, NORMALIZATION AND LOG2 FOLD CHANGE ####
# load pcr data
pcrRaw <- read.csv(paste0(file_path, pcrFile), row.names = 'Gene')
pcrRaw <- as.matrix(pcrRaw)
for(i in 1:nrow(pcrRaw)){
  if(length(which(pcrRaw[i,]==0))){
    pcrRaw[i,][which(pcrRaw[i,]==0)] <- mean(pcrRaw[i,][which(pcrRaw[i,]!=0)])
  }
}

# normalization
expData <- 2^(-pcrRaw)
sizeFactor <- calcNormFactors(expData, method = "TMM")
normData <- expData
for(i in 1:ncol(normData)) {
  normData[,i] <- expData[,i]/(sum(expData[,i])*sizeFactor[i])
}

# log2 folc change relative to sample mean
fcData <-t(apply(normData, 1, function(x) x/mean(x)))
log2fcData <-log2(fcData)


#### GENE EXPRESSION HEATMAP ####
# load gene expression data
heatData<-log2fcData
colnames(heatData)<-sapply(strsplit(colnames(heatData),'\\.'),function(x)x[2])

# column annotations
if(site){
  heatSite<-siteData
  names(heatSite)<-colnames(heatData)
  ha_col<-HeatmapAnnotation(Site=heatSite,
                            col=ann_colors[1],
                            annotation_name_side="left")
}else{
  ha_col<-NULL
}

# symetrical fold change color ramp (based on distribution of data, instead of min and max (which are often extreme outliers)
lowerBreak<-qnorm(0.01, mean(heatData), sd(heatData))
upperBreak<-qnorm(0.99, mean(heatData), sd(heatData))
color_fun<- colorRamp2(c(lowerBreak, mean(heatData), upperBreak), c("green4","black", "red3"))

# plot heatmap
ht1<-Heatmap(as.matrix(heatData),name = "log2 Fold Change",
             top_annotation = ha_col,
             show_row_dend=F,
             clustering_distance_columns = "manhattan",
             clustering_distance_rows = "manhattan",
             column_title = paste0(species, ' Gene Expression'),
             column_names_gp = gpar(fontsize = 10),
             row_names_gp = gpar(fontsize = 8),
             col = color_fun
)
ht1

#### GENE EXPRESSION PCA ####
# load log2FC data
pcaData <- data.frame(t(log2fcData), check.names = FALSE)
rownames(pcaData)<-sapply(strsplit(rownames(pcaData),'\\.'),function(x)x[2])

# PCA mod
pcaMod<-prcomp(pcaData, center=TRUE,scale=TRUE)

# Top 5 gene loading
PCA<-pcaMod
top5_pc1<-names(PCA$rotation[order(abs(PCA$rotation[,1]), decreasing=TRUE),1])[1:5]
top5_pc2<-names(PCA$rotation[order(abs(PCA$rotation[,2]), decreasing=TRUE),2])[1:5]
top5<-unique(c(top5_pc1, top5_pc2))

PCA2<-PCA
PCA2$rotation<-PCA2$rotation[top5,]
PCA2$center<-PCA2$center[top5]
PCA2$scale<-PCA2$scale[top5]

# Site option
if(site){
  pcaData$siteData<-siteData
}


p<-autoplot(PCA2, data=pcaData,label=F, 
            # optional arguments for colony of origin info
            colour=ifelse(site, 'siteData', 'red'),
            frame=ifelse(site, TRUE, FALSE), 
            
            # remaining common arguments
            shape = 16, size = 4,
            loadings = TRUE, loadings.label = TRUE,
            loadings.colour = 'darkgrey',
            loadings.label.colour = 'darkgrey',
            loadings.label.size = 4,
            main=paste0('PCA of ',species,' gene expressions')) + 
  geom_text(aes(label = rownames(pcaData)), size = 4,
            nudge_y = rep(0.025,rownames(pcaData)[1])) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black",fill=NA)) + 
  geom_hline(aes(yintercept=0), colour="black", linetype="dashed") +  
  geom_vline(aes(xintercept=0), colour="black", linetype="dashed")
p


#write.csv(cbind(top5_pc1,top5_pc2),'top_PCA of BLGU gene expressions.csv')


#### DISCRIMINANT ANALYSIS (COLONY OF ORIGIN, BLGU ONLY) ####

if(site){
  
  #Chem Data model
  modData<-mdl_filter(correctedData, mdl, looseFilter)
  modData<-data.frame(t(log10(modData)), check.names = FALSE)
  modData$colony<-siteData[-missingChem]
  #modData$colony<-batchData[-missingChem]
  plsMod<-plsda(select(modData, !colony), as.factor(modData$colony),
                center=T,
                scale=T,
                cv=1,
                #cv = list("rand", 4, 20),
                ncomp.selcrit = "min")
  
  plot(plsMod, labels=TRUE)
  summary(plsMod)
  
  #ge Data model
  modData<-data.frame(t(log2fcData),check.names = FALSE)
  modData$colony<-siteData
  #modData$colony<-batchData
  plsMod<-plsda(select(modData, !colony), as.factor(modData$colony),
                center=T,
                scale=T,
                cv=1,
                #cv = list("rand", 4, 20),
                ncomp.selcrit = "min")
  
  plot(plsMod, labels=TRUE)
  summary(plsMod)
}



#### PLS REGRESSION (separate PLS models. DO NOT USE) ####

# chem data for pls
plsChemData<-mdl_filter(correctedData, mdl, 0.75)
plsChemData<-data.frame(t(log10(plsChemData)), check.names = FALSE)

# gene expression data for pls
plsGeneData<-t(log2fcData[,-missingChem])

# results storage vectors
nGenes<-vector()
nComps<-vector()
xVar<-vector()
yVar<-vector()
sMean<-vector()
sMed<-vector()
sMin<-vector()
sMax<-vector()

#pls model-all genes
plsMod<-pls(plsGeneData, plsChemData, center=TRUE, scale=TRUE, cv=1)
nGenes<-c(nGenes,ncol(plsGeneData))
nComps<-c(nComps,plsMod$ncomp.selected)

#plot top slopes
par(mfrow=c(3,3), mar=c(4,3,1.5,1))
for(i in order(plsMod$cvres$slope[,3],decreasing=TRUE)[1:9]){
  plot(
    plsMod$cvres$y.pred[,3,i]~plsMod$cvres$y.ref[,i],
    xlab=NA,
    ylab=NA,
    main=dimnames(plsMod$calres$y.pred)[[3]][i],
    col='red', pch=16, cex=1.8
  )
  text(plsMod$cvres$y.ref[,i], plsMod$cvres$y.pred[,nComps,i], labels=rownames(plsChemData))
  tempLM<-lm(plsMod$cvres$y.pred[,nComps,i]~plsMod$cvres$y.ref[,i])
  abline(tempLM, col="blue", lwd=2)
  abline(0,1, lty=2)
  legend('bottomright', legend = round(plsMod$cvres$slope[i,nComps],2))
}
plsMod.2 <- plsMod


rmseMean<-mean(plsMod.2$cvres$rmse[,nComps[1]])
rmseMed<-median(plsMod.2$cvres$rmse[,nComps[1]])
rmseMin<-min(plsMod.2$cvres$rmse[,nComps[1]])
rmseMax<-max(plsMod.2$cvres$rmse[,nComps[1]])

xVar<-plsMod.2$calres$xdecomp$cumexpvar[nComps[1]]
yVar<-plsMod.2$calres$ydecomp$cumexpvar[nComps[1]]


sMean<-c(sMean, mean(plsMod.2$cvres$slope[,1]))
sMed<-c(sMed, median(plsMod.2$cvres$slope[,1]))
sMin<-c(sMin, min(plsMod.2$cvres$slope[,1]))
sMax<-c(sMax, max(plsMod.2$cvres$slope[,1]))


#pls model-specific VIP gene for each chemical
plsVIPscores<-vipscores(plsMod.2)
hasVIPs<-names(which(apply(plsVIPscores, 2, function(x){any(x>1)})))
#unique Mods based on VIPs per chem
uniqueMods<-list()
uniqueVIPgene<-list()

for(i in hasVIPs){
  uniqueVIPgene[[i]]<-plsGeneData[,plsVIPscores[,i]>1]
  uniqueMods[[i]]<-pls(as.matrix(uniqueVIPgene[[i]]), as.matrix(plsChemData[,i]), center=TRUE, scale=TRUE, cv=1)
}
plotOrder<-order(sapply(uniqueMods, function(x){x$cvres$slope[,x$ncomp.selected]}), decreasing = TRUE)


nGenes<-c(nGenes,length(hasVIPs))
nComps<-c(nComps,1)

rmseMean<-c(rmseMean, mean(sapply(uniqueMods, function(x) x$cvres$rmse[,nComps[1]])))
rmseMed<-c(rmseMed, median(sapply(uniqueMods, function(x) x$cvres$rmse[,nComps[1]])))
rmseMin<-c(rmseMin, min(sapply(uniqueMods, function(x) x$cvres$rmse[,nComps[1]])))
rmseMax<-c(rmseMax, max(sapply(uniqueMods, function(x) x$cvres$rmse[,nComps[1]])))

tempVar<-median(sapply(uniqueMods, function(x) x$calres$xdecomp$cumexpvar[nComps[1]]))
xVar<-c(xVar, tempVar)

tempVar<-median(sapply(uniqueMods, function(x) x$calres$ydecomp$cumexpvar[nComps[1]]))
yVar<-c(yVar, tempVar)

sMean<-c(sMean, mean(sapply(uniqueMods, function(x) x$cvres$slope[,1])))
sMed<-c(sMed, median(sapply(uniqueMods, function(x) x$cvres$slope[,1])))
sMin<-c(sMin, min(sapply(uniqueMods, function(x) x$cvres$slope[,1])))
sMax<-c(sMax, max(sapply(uniqueMods, function(x) x$cvres$slope[,1])))


par(mfrow=c(4,3), mar=c(4,3,1.5,1))
for(i in plotOrder){
  u_nComps<-uniqueMods[[i]]$ncomp.selected
  dataRange<-c(
    min(c(uniqueMods[[i]]$cvres$y.pred[,u_nComps,1],uniqueMods[[i]]$cvres$y.ref[,1]))-0.2,
    max(c(uniqueMods[[i]]$cvres$y.pred[,u_nComps,1],uniqueMods[[i]]$cvres$y.ref[,1]))+0.2
  )
  plot(
    uniqueMods[[i]]$cvres$y.pred[,u_nComps,1]~uniqueMods[[i]]$cvres$y.ref[,1],
    xlab=NA,
    ylab=NA,
    main=paste0(colnames(plsChemData)[i]),
    ylim=dataRange,
    xlim=dataRange,
    col='red', pch=16, cex=1.8
  )
  print(colnames(plsChemData)[i])
  #slope <- round(uniqueMods[[i]]$cvres$slope[,uniqueMods[[i]]$ncomp.selected],2)
  #mtext(side=1, line=2, cex=0.7, paste0('slope = ', slope, adj=0))
  #r2<-format(round(uniqueMods[[i]]$cvres$r2[,uniqueMods[[i]]$ncomp.selected],2))
  #mtext(side=1, line=2, cex=0.7, bquote(r^2 == .(r2)), adj=1)
  mtext(paste0(
    "slope=", round(uniqueMods[[i]]$cvres$slope[,u_nComps],2),
    " rmse=", round(uniqueMods[[i]]$cvres$rmse[,u_nComps],2)
  ),
  side=1, adj=1, line=2, cex=0.8)
  #mtext(side=1, line=2, cex=0.7, bquote(r^2 == .(r2)), adj=1)
  #text(uniqueMods[[i]]$cvres$y.ref[,1], uniqueMods[[i]]$cvres$y.pred[,u_nComps,1], labels=rownames(plsChemData))
  tempLM<-lm(uniqueMods[[i]]$cvres$y.pred[,u_nComps,1]~uniqueMods[[i]]$cvres$y.ref[,1])
  abline(tempLM, col="blue", lwd=2)
  abline(0,1, lty=2)
}


###top genes
VIPchem <- plsVIPscores[,hasVIPs]
topgene <- matrix(nrow=10, ncol=ncol(VIPchem))
colnames(topgene) <- colnames(VIPchem)
write.csv(apply(VIPchem, 2, function(x) names(x[order(x)[length(x):1]][1:10])),
          ' top gene list-BLGU.csv')

### summary of variance
modRes2<-data.frame(
  model=1:2,
  type=c("all genes","VIP gene >1"),
  nGenes=nGenes,
  nComps=nComps,
  gene_Var=xVar,
  chem_Var=yVar,
  rmseMean=rmseMean,
  rmseMed=rmseMed,
  rmseMin=rmseMin,
  rmseMax=rmseMax,
  slopeMean=sMean,
  slopeMed=sMed,
  slopeMin=sMin,
  slopeMax=sMax,
  stringsAsFactors=FALSE
)
modRes2


#### PLS REGRESSION (sinlge PLS model. USE THIS) ####

# chem data for pls
plsChemData<-mdl_filter(correctedData, mdl, stringentFilter)
plsChemData<-data.frame(t(log10(plsChemData)), check.names = FALSE)

# gene expression data for pls
plsGeneData<-t(log2fcData[,-missingChem])

# results storage vectors
nGenes<-vector()
nComps<-vector()
xVar<-vector()
yVar<-vector()
sMean<-vector()
sMed<-vector()
sMin<-vector()
sMax<-vector()

#pls model-all genes
plsMod<-pls(plsGeneData, plsChemData, center=TRUE, scale=TRUE, cv=1)
nGenes<-c(nGenes,ncol(plsGeneData))
nComps<-c(nComps,plsMod$ncomp.selected)

#plot top slopes
par(mfrow=c(3,3), mar=c(4,3,1.5,1))
for(i in order(plsMod$cvres$slope[,nComps[1]],decreasing=TRUE)[1:9]){
  dataRange<-c(
    min(c(plsMod$cvres$y.pred[,nComps[1],i],plsMod$cvres$y.ref[,i]))-0.2,
    max(c(plsMod$cvres$y.pred[,nComps[1],i],plsMod$cvres$y.ref[,i]))+0.2
  )
  
  plot(
    plsMod$cvres$y.pred[,nComps[1],i]~plsMod$cvres$y.ref[,i],
    xlab=NA,
    ylab=NA,
    ylim=dataRange,
    xlim=dataRange,
    main=dimnames(plsMod$calres$y.pred)[[3]][i],
    col='red', pch=16, cex=1.8
  )
  text(plsMod$cvres$y.ref[,i], plsMod$cvres$y.pred[,nComps[1],i], labels=rownames(plsChemData))
  tempLM<-lm(plsMod$cvres$y.pred[,nComps[1],i]~plsMod$cvres$y.ref[,i])
  abline(tempLM, col="blue", lwd=2)
  abline(0,1, lty=2)
  mtext(paste0(
    "slope=", round(plsMod$cvres$slope[i,nComps[1]],2),
    " rmse=", round(plsMod$cvres$rmse[i,nComps[1]],2)
  ),
  side=1, adj=1, line=2, cex=0.8)
  
}


#pls with VIP genes only
plsVIPscores<-vipscores(plsMod)
hasVIPs<-names(which(apply(plsVIPscores, 2, function(x){any(x>1)})))

#pls with VIP>1
exclrows <- apply(plsVIPscores,2,function(x){all(x<1)})
exclcols <- apply(plsVIPscores,1,function(x){all(x<1)})
plsV2 = pls(plsGeneData, plsChemData,
            center=TRUE, scale=TRUE, cv = 1,
            exclrows = exclrows, exclcols = exclcols
)

# plot results
v2Comps<-plsV2$ncomp.selected

par(mfrow=c(3,3), mar=c(4,3,1.5,1))
for(i in order(plsV2$cvres$slope[,v2Comps[1]],decreasing=TRUE)[1:length(hasVIPs)]){
  dataRange<-c(
    min(c(plsV2$cvres$y.pred[,v2Comps[1],i],plsV2$cvres$y.ref[,i]))-0.2,
    max(c(plsV2$cvres$y.pred[,v2Comps[1],i],plsV2$cvres$y.ref[,i]))+0.2
  )
  
  plot(
    plsV2$cvres$y.pred[,v2Comps[1],i]~plsV2$cvres$y.ref[,i],
    xlab=NA,
    ylab=NA,
    ylim=dataRange,
    xlim=dataRange,
    main=dimnames(plsV2$calres$y.pred)[[3]][i],
    col='red', pch=16, cex=1.8
  )
  text(plsV2$cvres$y.ref[,i], plsV2$cvres$y.pred[,v2Comps[1],i], labels=rownames(plsChemData))
  tempLM<-lm(plsV2$cvres$y.pred[,v2Comps[1],i]~plsV2$cvres$y.ref[,i])
  abline(tempLM, col="blue", lwd=2)
  abline(0,1, lty=2)
  mtext(paste0(
    "slope=", round(plsV2$cvres$slope[i,v2Comps[1]],2),
    " rmse=", round(plsV2$cvres$rmse[i,v2Comps[1]],2)
  ),
  side=1, adj=1, line=2, cex=0.8)
  
}








####################################################################################
### sum by group analysis
####################################################################################
merge.data<-cbind(rownames(correctedData),correctedData)
colnames(merge.data)[1]<-colnames(group)[1]
merge.data<-merge(group,merge.data,by=colnames(group)[1])
TE<-which(merge.data[,3]=='Trace Elements') # here we dont want to sum trace elements by groups
PAC_sum<-apply(merge.data[-TE,-c(1:3)],2,function(x) {tapply(x,merge.data[-TE,2],sum)})
PAC_sum<-PAC_sum[-which(rownames(PAC_sum)=='Trace Elements'),]
TE<-merge.data[TE,-(2:3)]
rownames(TE)<-TE[,1]
TE<-TE[,-1]
sum.data<-rbind(TE,PAC_sum)

heatData<-sum.data
heatData<-log10(heatData)
colnames(heatData)<-sapply(strsplit(colnames(heatData),'\\.'),function(x)x[2])

heatGroup<-c(rep('Trace Elements',5),rownames(heatData)[-c(1:5)])

# row annotations
Group<-factor(heatGroup, levels=unique(heatGroup))

ha_row<-rowAnnotation(Group=Group,
                      col=ann_colors[2],
                      annotation_name_side="top")

# column annotations (if site=TRUE, or batch=TRUE)
if(site){
  heatSite<-siteData[-missingChem]
  names(heatSite)<-colnames(heatData)
  ha_col<-HeatmapAnnotation(Site=heatSite,
                            col=ann_colors[1],
                            annotation_name_side="left")
}else{
  ha_col<-NULL
}

if(batch){
  heatBatch<-batchData[-missingChem]
  names(heatBatch)<-colnames(heatData)
  ha_col<-HeatmapAnnotation(
    Batch=heatBatch,
    col=ann_colors[4],
    annotation_name_side="left")
}

if(site&batch){
  ha_col<-HeatmapAnnotation(
    Batch=heatBatch,
    Site=heatSite,
    col=ann_colors[c(1,4)],
    annotation_name_side="left")
}



# concentration colour palette
color_fun<- colorRamp2(c(min(heatData),0,max(heatData)), c("cornflowerblue", "#EEEEEE", "red"))

#plot heatmap
ht1<-Heatmap(as.matrix(heatData),name = "log10 (ng/g)",
             top_annotation = ha_col,
             right_annotation = ha_row,
             show_row_dend=F,
             clustering_distance_columns = "manhattan",
             column_title = paste0(species, ' Chemical Residues'),
             cluster_rows = F,
             column_names_gp = gpar(fontsize = 10),
             row_names_gp = gpar(fontsize = 8),
             col = color_fun)
ht1

#### CHEMICAL PCA ####
# filter and transcform chemical data
pcaData<-sum.data
pcaData<-t(log10(pcaData))
pcaData<-data.frame(pcaData, check.names =FALSE)
rownames(pcaData)<-sapply(strsplit(rownames(pcaData),'\\.'),function(x)x[2])

# pca model
pcaMod<-prcomp(pcaData,center=TRUE)

# top loadings
top5_pc1<-names(pcaMod$rotation[order(abs(pcaMod$rotation[,1]), decreasing=TRUE),1])[1:5]
top5_pc2<-names(pcaMod$rotation[order(abs(pcaMod$rotation[,2]), decreasing=TRUE),2])[1:5]
top5<-unique(c(top5_pc1, top5_pc2))

pcaMod2<-pcaMod
pcaMod2$rotation<-pcaMod2$rotation[top5,]
pcaMod2$center<-pcaMod2$center[top5]
pcaMod2$scale<-pcaMod2$scale[top5]


# pca color and frame options
pcaCol<-ifelse(site|batch, "legend", "red")
pcaFrame<-ifelse(site, TRUE, ifelse(batch, TRUE , FALSE))  

if(site){
  pcaData$legend<-siteData[-missingChem]
}

if(batch){
  pcaData$legend<-batchData[-missingChem]
}

# plot pca
p<-autoplot(pcaMod2, 
            # optional arguments for colony of origin info
            colour=pcaCol,
            frame=pcaFrame,
            # remaining common arguments
            shape = 16, size = 4,
            data=pcaData,label=F, 
            loadings = TRUE, loadings.label = TRUE,
            loadings.colour = 'darkgrey',
            loadings.label.colour = 'darkgrey',
            loadings.label.size = 4,
            main=paste0('PCA of ',species,' chemical burdens')) + 
  geom_text(aes(label = rownames(pcaData)), size = 4,
            nudge_y = rep(0.025,dim(pcaData)[1])) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black",fill=NA))
p



#### PLS REGRESSION (sinlge PLS model. USE THIS) ####

# chem data for pls
plsChemData<-sum.data
plsChemData<-data.frame(t(log10(plsChemData)), check.names = FALSE)

# gene expression data for pls
plsGeneData<-t(log2fcData[,-missingChem])

# results storage vectors
nGenes<-vector()
nComps<-vector()
xVar<-vector()
yVar<-vector()
sMean<-vector()
sMed<-vector()
sMin<-vector()
sMax<-vector()

#pls model-all genes
plsMod<-pls(plsGeneData, plsChemData, center=TRUE, scale=TRUE, cv=1)
nGenes<-c(nGenes,ncol(plsGeneData))
nComps<-c(nComps,plsMod$ncomp.selected)

#plot top slopes
par(mfrow=c(3,3), mar=c(4,3,1.5,1))
for(i in order(plsMod$cvres$slope[,nComps[1]],decreasing=TRUE)[1:9]){
  dataRange<-c(
    min(c(plsMod$cvres$y.pred[,nComps[1],i],plsMod$cvres$y.ref[,i]))-0.2,
    max(c(plsMod$cvres$y.pred[,nComps[1],i],plsMod$cvres$y.ref[,i]))+0.2
  )
  
  plot(
    plsMod$cvres$y.pred[,nComps[1],i]~plsMod$cvres$y.ref[,i],
    xlab=NA,
    ylab=NA,
    ylim=dataRange,
    xlim=dataRange,
    main=dimnames(plsMod$calres$y.pred)[[3]][i],
    col='red', pch=16, cex=1.8
  )
  text(plsMod$cvres$y.ref[,i], plsMod$cvres$y.pred[,nComps[1],i], labels=rownames(plsChemData))
  tempLM<-lm(plsMod$cvres$y.pred[,nComps[1],i]~plsMod$cvres$y.ref[,i])
  abline(tempLM, col="blue", lwd=2)
  abline(0,1, lty=2)
  mtext(paste0(
    "slope=", round(plsMod$cvres$slope[i,nComps[1]],2),
    " rmse=", round(plsMod$cvres$rmse[i,nComps[1]],2)
  ),
  side=1, adj=1, line=2, cex=0.8)
  
}


#pls with VIP genes only
plsVIPscores<-vipscores(plsMod)
hasVIPs<-names(which(apply(plsVIPscores, 2, function(x){any(x>1)})))

#pls with VIP>1
exclrows <- apply(plsVIPscores,2,function(x){all(x<1)})
exclcols <- apply(plsVIPscores,1,function(x){all(x<1)})
plsV2 = pls(plsGeneData, plsChemData,
            center=TRUE, scale=TRUE, cv = 1,
            exclrows = exclrows, exclcols = exclcols
)

# plot results
v2Comps<-plsV2$ncomp.selected

par(mfrow=c(3,3), mar=c(4,3,1.5,1))
for(i in order(plsV2$cvres$slope[,v2Comps[1]],decreasing=TRUE)[1:length(hasVIPs)]){
  dataRange<-c(
    min(c(plsV2$cvres$y.pred[,v2Comps[1],i],plsV2$cvres$y.ref[,i]))-0.2,
    max(c(plsV2$cvres$y.pred[,v2Comps[1],i],plsV2$cvres$y.ref[,i]))+0.2
  )
  
  plot(
    plsV2$cvres$y.pred[,v2Comps[1],i]~plsV2$cvres$y.ref[,i],
    xlab=NA,
    ylab=NA,
    ylim=dataRange,
    xlim=dataRange,
    main=dimnames(plsV2$calres$y.pred)[[3]][i],
    col='red', pch=16, cex=1.8
  )
  text(plsV2$cvres$y.ref[,i], plsV2$cvres$y.pred[,v2Comps[1],i], labels=rownames(plsChemData))
  tempLM<-lm(plsV2$cvres$y.pred[,v2Comps[1],i]~plsV2$cvres$y.ref[,i])
  abline(tempLM, col="blue", lwd=2)
  abline(0,1, lty=2)
  mtext(paste0(
    "slope=", round(plsV2$cvres$slope[i,v2Comps[1]],2),
    " rmse=", round(plsV2$cvres$rmse[i,v2Comps[1]],2)
  ),
  side=1, adj=1, line=2, cex=0.8)
  
}
