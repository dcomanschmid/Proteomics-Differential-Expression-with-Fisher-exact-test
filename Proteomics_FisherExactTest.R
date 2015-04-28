#####################################################################################################                                                                                                
# Proteomics data analysis: Fisher exact test for differential expression                           #
#	Input:                                                                                      #
#		- table (TAB delimited) as .txt, .csv, .xlsx file                                   #
#		- rows=spectral counts                                                              #
#		- columns=conditions                                                                #
#	Analysis:                                                                                   #
#		- define (all possible) pairwise comparisons                                        #
#		- filter proteins with low/high spectral counts                                     #
#		- Fisher exact test                                                                 #
#		- p-value and multiple testing correction                                           #
# 	Output:               								            #
#		- Excel workbook with one sheet for each pairwise comparison                        #
#         	- sheet content: ProtID,Cond1,Cond2,log2FC,p-value, p-adjust (annotation, etc.)     #
#####################################################################################################

library("stats")
library("reshape")
library("XLConnect")
library ("gplots")

    # define the working directory 

datpath = "/workingDir/ctrl1x/"

    # read in the input data available as one folder per sample (reference or control, treatment); each folder contains one file per technical replicate
    # file formats: .txt, .csv etc.

filelist = list.files(datpath, pattern = "*.csv")

    # load the columns with the unique spectral counts [e.g. c(1,3)] or all spectral counts [e.g.c(1,2)]  

datalist = lapply(filelist, function(x)read.csv(file.path(paste(datpath,x,sep="")), header=F,sep=";")[,c(1,2)])

    # rename the columns in each file     

names(datalist) <- filelist
for (d in names(datalist)){
  colnames(datalist[[d]]) <- c("protID","SC")
}

summary(datalist)


    # ALTERNATIVE: read in the input data available as one file for all samples
    # file formats: .txt, .csv etc.

# sp.counts <- read.table(file.path(paste(datpath,"xxxxx.txt",sep="")),row.names=1,header=TRUE)

    # load protein annotation 

# prot.anno <- read.csv(file.path(paste(datpath,"ProtAnno_Isabel.csv",sep="")),sep=";",header=TRUE)
# colnames(prot.anno) <- c("Row.names","Prot.anno")


    # make data matrices for reference (control) and treatment samples by summing the spectral counts from technical replicates
    # see TO DO above
    
ctrl1x <- Reduce(function(x,y) merge(x,y, all=T,by.x='protID',by.y='protID'),datalist, accumulate=F)
colnames(ctrl1x) <- c("protID",filelist)
ctrl1x$sum <- rowSums(ctrl1x[,2:ncol(ctrl1x)],na.rm=T)


ctrl10x <- Reduce(function(x,y) merge(x,y, all=T,by.x='protID',by.y='protID'),datalist, accumulate=F)
colnames(ctrl10x) <- c("protID",filelist)
ctrl10x$sum <- rowSums(ctrl10x[,2:ncol(ctrl10x)],na.rm=T)

mt1x <- Reduce(function(x,y) merge(x,y, all=T,by.x='protID',by.y='protID'),datalist, accumulate=F)
colnames(mt1x) <- c("protID",filelist)
mt1x$sum <- rowSums(mt1x[,2:ncol(mt1x)],na.rm=T)

mt10x <- Reduce(function(x,y) merge(x,y, all=T,by.x='protID',by.y='protID'),datalist, accumulate=F)
colnames(mt10x) <- c("protID",filelist)
mt10x$sum <- rowSums(mt10x[,2:ncol(mt10x)],na.rm=T)

    # collect data in one list and make a data frame with rows=proteins and columns=samples (technical replicates summed)

allsamp.sum <- list(ctrl1x[,c("protID","sum")],ctrl10x[,c("protID","sum")],mt1x[,c("protID","sum")],mt10x[,c("protID","sum")])
sc.data <- Reduce(function(x,y) merge(x,y, all=T,by.x='protID',by.y='protID'),allsamp.sum, accumulate=F)
sc.data[is.na(sc.data)] <- 0
colnames(sc.data) <- c("Row.names","ctrl1x","ctrl10x","mt1x","mt10x")


    # make a data frame for each (possible) pairwise comparison
    # filter proteins with spectral  counts <2 or > 5000 
    # store the filtered pairwise data frames in a list

pw <- combn(ncol(sc.data),m=2)
dim(pw)
pw
pw.l <- list()
for (p in 1:ncol(pw)){
  pw.df <- sc.data[,pw[,p]]
  pw.df.filt <- pw.df[rowSums(pw.df) >= 2 & rowSums(pw.df) <= 5000,]
  pw.l[[p]] <- pw.df.filt
}

names.pw <- combn(colnames(sc.data),m=2)
names.pwfull <- character()
for (n in 1:ncol(names.pw)){
  names.pwfull <- c(names.pwfull,paste(names.pw[1,n],names.pw[2,n],sep="_"))
}
names.pwfull

names(pw.l) <- names.pwfull
summary(pw.l)


    # plot the number of proteins after filtering

filt.p <- numeric()
for (n in names(pw.l)){
  filt.p[n] <- nrow(pw.l[[n]]) 
}

par(mar=c(2,11,4,2))
barplot(filt.p,col="aliceblue",main="#Prot. with SC >=2 & <= 5000",horiz=T,las=1)


    # make the contingency table (for Fisher exact test) for each pariwise comparison


ctpw.l <- list()
for (l in names(pw.l)){
  ct.l <- list()
  for (p in seq(along=pw.l[[l]][,1])){
    a <- pw.l[[l]][p,1]
    b <- sum(pw.l[[l]][,1])- pw.l[[l]][p,1]
    cc <- pw.l[[l]][p,2]
    d <- sum(pw.l[[l]][,2])-pw.l[[l]][p,2]
    ct.l[[p]] <- cbind(a,b,cc,d)
    names(ct.l)[p] <- row.names(pw.l[[l]])[p]
  }
  ctpw.l[[l]] <- ct.l
}

summary(ctpw.l)


    # apply the Fisher exact test for each protein in each pariwise comparison
    # extract the p-value 
    # apply multiple testing correction (FDR or other)

pvalpw.l <- list()
for (q in names(ctpw.l)){
  pval.l <- list()
  for (i in names(ctpw.l[[q]])){
    fisherpval <- apply(ctpw.l[[q]][[i]],1, function(x) fisher.test(matrix(x,nr=2),simulate.p.value=TRUE)$p.value)
    fisherpval_df <- as.data.frame(fisherpval)
#   row.names(fisherpval_df) <- names(ct.l)
    pval.l[[i]] <- fisherpval_df
  }
  pval.all <-do.call('rbind',pval.l)
  pval.all$padj <- p.adjust(pval.all$fisherpval,"fdr", length(pval.all$fisherpval))
  pvalpw.l[[q]] <- pval.all
}

summary(pvalpw.l)

  # add log2 FC column 

diffE <- list()
for (f in names(pw.l)){
  diffE[[f]] <- merge(pw.l[[f]],pvalpw.l[[f]],by="row.names")
  diffE[[f]]$FC <- diffE[[f]][,3] / diffE[[f]][,2]
  diffE[[f]]$log2FC <- log2(diffE[[f]]$FC)
}

summary(diffE)


  # add protein annotation
  
diffE.anno <- list()
for (d in names(diffE)){
  diffE.anno[[d]] <- merge(diffE[[d]],prot.anno,by="Row.names")
}

summary(diffE.anno)


  # write results in an Excel workbook (one sheet for each pairwsie comparison)
  # sheet content: ProtID,Cond1,Cond2,log2FC,p-values and p-adjust, annotation       
  
wb.res <- loadWorkbook(file.path(paste(datpath,"xxxxx_SumTechRepl_FisherDiffEProt.xlsx",sep="")), create = TRUE)
createSheet(wb.res, name=names(diffE.anno))
writeWorksheet(wb.res, diffE, names(diffE.anno),header=T)
saveWorkbook(wb.res)

    
# TO DO: (automatic read in data from different directories and sum the spectral counts from technical replicates): 
#     
# 	fd <- c("ctrl1x","ctrl10x","mt1x","mt10x")
# 	
#       sc <- list()
#       for (f in 1:length(fd)){
# 	     
#            datpath = "/.../f/"
# 	     filelist = list.files(datpath, pattern = "*.csv")
# 	     datalist = lapply(filelist, function(x)read.csv(file.path(paste(datpath,x,sep="")), header=F,sep=";")[,c(1,2)])
#           
#            names(datalist) <- filelist
# 	     for (d in names(datalist)){
#                 colnames(datalist[[d]]) <- c("protID","SC")
#            }
#            
#            sc[[f]] <- Reduce(function(x,y) merge(x,y, all=T,by.x='protID',by.y='protID'),datalist, accumulate=F)
#            colnames(sc)[[f]] <- c("protID",filelist)
#            sc[[f]]$sum <- rowSums(sc[[f]][,2:ncol(sc[[f]])],na.rm=T)
#	     sc[[f]] <- sc[[f]][,c("protID","sum")]
#        }
# sc.data <- Reduce(function(x,y) merge(x,y, all=T,by.x='protID',by.y='protID'),sc, accumulate=F)
# sc.data[is.na(sc.data)] <- 0
# colnames(sc.data) <- c("Row.names","ctrl1x","ctrl10x","mt1x","mt10x")
