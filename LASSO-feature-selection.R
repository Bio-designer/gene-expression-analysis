#############################
######Feature selection by LASSO in gene expression data
######@Minstein, 2017-07-07
#############################

##load the glmnet package
library(glmnet)

input.dir="path of your input datasets" ##please edit this before test the program


##function copied from MP.Library.R (http://portals.broadinstitute.org/cgi-bin/cancer/publications/view/161).
MP.Gct2Frame <- function(filename = "NULL") {
  #
  # Reads a gene expression dataset in GCT format and converts it into an R data frame
  #
  # The Broad Institute
  # SOFTWARE COPYRIGHT NOTICE AGREEMENT
  # This software and its documentation are copyright 2003 by the
  # Broad Institute/Massachusetts Institute of Technology.
  # All rights are reserved.
  #
  # This software is supplied without any warranty or guaranteed support
  # whatsoever. Neither the Broad Institute nor MIT can be responsible for
  # its use, misuse, or functionality.

  ds <- read.delim(filename, header=T, sep="\t", skip=2, row.names=1, blank.lines.skip=T, comment.char="", as.is=T, na.strings = "")
  descs <- ds[,1]
  ds <- ds[-1]
  row.names <- row.names(ds)
  names <- names(ds)
  return(list(ds = ds, row.names = row.names, descs = descs, names = names))
}

##function copied from MP.Library.R (http://portals.broadinstitute.org/cgi-bin/cancer/publications/view/161).
MP.ReadClsFile <- function(file = "NULL") {
  #
  # Reads a class vector CLS file and defines phenotype and class labels vectors (numeric and character) for the samples in a gene expression file (RES or GCT format)
  #
  # The Broad Institute
  # SOFTWARE COPYRIGHT NOTICE AGREEMENT
  # This software and its documentation are copyright 2003 by the
  # Broad Institute/Massachusetts Institute of Technology.
  # All rights are reserved.
  #
  # This software is supplied without any warranty or guaranteed support
  # whatsoever. Neither the Broad Institute nor MIT can be responsible for
  # its use, misuse, or functionality.

  cls.cont <- readLines(file)
  num.lines <- length(cls.cont)
  class.list <- unlist(strsplit(cls.cont[[3]], " "))
  s <- length(class.list)
  t <- table(class.list)
  l <- length(t)
  phen <- vector(length=l, mode="character")
  class.v <- vector(length=s, mode="numeric")

  current.label <- class.list[1]
  current.number <- 1
  class.v[1] <- current.number
  phen[1] <- current.label
  phen.count <- 1

  if (length(class.list) > 1) {
    for (i in 2:s) {
      if (class.list[i] == current.label) {
        class.v[i] <- current.number
      } else {
        phen.count <- phen.count + 1
        current.number <- current.number + 1
        current.label <- class.list[i]
        phen[phen.count] <- current.label
        class.v[i] <- current.number
      }
    }
  }
  return(list(phen = phen, class.v = class.v, class.list = class.list))
}

############################



# Read dataset using function MP.Gct2Frame()
# Samples from Ross et al 2003 (PMID: 12730115) and Ross et al 2004 (PMID: 15226186)
gct.file=paste(input.dir,"Leukemia2.all.gct",sep="")
dataset<-MP.Gct2Frame (gct.file)
##raw matrix
rm <- data.matrix(dataset$ds)
gs.names <- dataset$row.names
gs.descs <- dataset$descs
sample.names <- dataset$names

#Read CLS file using function MP.ReadClsFile()
cls.file<-paste(input.dir,"Leukemia2.all.cls",sep="")
CLS <- MP.ReadClsFile(file=cls.file)
class.labels <- CLS$class.v
##phenotype
class.phen <- CLS$phen
class.list <- CLS$class.list


##Divide into three major classes
#ALL-T
#ALL.T. T cell lineage leukemia
ALLT.pos<-grep("ALL.T...",class.list,fixed=TRUE)
#ALLT.pos<-which(class.list%in%c("ALL.T.....D1",""))
ALLT.count<-length(ALLT.pos)
#have a check
#class.list[ALLT.pos]

##AML
AML.pos<-grep("AML.",class.list,fixed=TRUE)
#AML.pos<-which(class.list%in%class.phen[6:10])
AML.count<-length(AML.pos)
#class.list[AML.pos]

#ALL-B
#B cell lineage leukemia
ALLB.pos<-seq(1,length(class.list),1)[-c(ALLT.pos,AML.pos)]
#ALLB.pos<-which(class.list%in%class.phen[c(1,2,3,5)])
ALLB.count<-length(ALLB.pos)
#class.list[ALLB.pos]


#####Generate trianing dataset and test dataset
getSample.pos<-function(pos1,pos2,pos3,size="NULL",seed){

  if(size[1]=="NULL"){
    size=c(10,10,10)
  }
  ##Is this useful for sampling?
  set.seed(seed)


  temp1<-sample(pos1,size[1])
  temp2<-sample(pos2,size[2])
  temp3<-sample(pos3,size[3])

  pos<-c(temp1,temp2,temp3)
  ##-1 for ALL-T, 0 for ALL-B
  classlabel<-c(rep(-1,size[1]),rep(0,size[2]),
                rep(1,size[3]))
  classname<-c(rep("ALL-T",size[1]),rep("ALL-B",size[2]),
               rep("AML",size[3]))

  return(list(pos=pos,classlabel=classlabel,classname=classname))
}

##function copied from MP.Library.R (http://portals.broadinstitute.org/cgi-bin/cancer/publications/view/161).

MP.NormalizeCols.Rank <- function(V) {

  cols <- length(V[1,])
  rows <- length(V[,1])
  for (j in 1:cols) {  # column rank normalization
    V[,j] <- rank(V[,j], ties.method = "average")
  }

  return(V)
}

##Preprocess the data
data.preprocess <- function(m) {

  ##Preprocess m
  # threshold, ceiling and shift

  ##Calculate the 0.05 quantile and 0.95 quantile first
  cols <- length(m[1,])
  quantile_up<-rep(0,cols)
  for (j in 1:cols) {
    quantile_up[j]<-quantile(m[,j],0.95)
  }
  quantile_low<-rep(0,cols)
  for (j in 1:cols) {
    quantile_low[j]<-quantile(m[,j],0.05)
  }
  #alternative method to calculate quantile
  #https://stats.stackexchange.com/questions/13399/calculating-the-95th-percentile-comparing-normal-distribution-r-quantile-and
  #sort(m[,1])[0.95*length(m[,1])]

  thres<-mean(quantile_low)## set by summary(m) first and to evaluate the numbers by which(m>ceil) and which(m<thres)
  ceil<-mean(quantile_up)

  m[m < thres] <- thres
  m[m > ceil] <- ceil

  ##this suits best now
  m <- MP.NormalizeCols.Rank(m)/length(m[,1])

  return(m)

}

ALLT.size<-ceil(length(ALLT.pos)*0.5)
ALLB.size<-ceil(length(ALLB.pos)*0.5)
AML.size<-ceil(length(AML.pos)*0.5)



##set the seed to 940309, size to c(32,32,32)
#train.table<-getTrian.pos(ALLT.pos,ALLB.pos,AML.pos,size = c(ALLT.size,ALLB.size,AML.size),seed=940309)
train.table<-getSample.pos(ALLT.pos,ALLB.pos,AML.pos,size = c(ALLT.size,ALLT.size,ALLT.size),seed=940309)
train.pos<-train.table$pos
train.classlabel<-train.table$classlabel
#train.classname<-class.list[train.pos]
train.classname<-train.table$classname


##Generate the training dataset
train<-rm[,train.pos]


##Preprocess the training dataset
train<-data.preprocess(train)


## cross-validation fit
cvfit<- cv.glmnet(t(train),train.classlabel)
#It includes the cross-validation curve (red dotted line), and upper and lower standard deviation curves along the λ sequence (error bars).
#Two selected λ’s are indicated by the vertical dotted lines.
plot(cvfit)
print(cvfit)


#We can view the selected λ’s and the corresponding coefficients.
#lambda.min is the value of λ that gives minimum mean cross-validated error.
cvfit$lambda.min

#The other λ saved is lambda.1se,
#which gives the most regularized model such that error is within one standard error of the minimum.
#To use that, we only need to replace lambda.min with lambda.1se above.
cvfit$lambda.1se

##use this to control the degree of freedom
coef_temp<-coef(cvfit, s = "lambda.1se")
str(coef_temp)


##locations for the selected genes of which the coefficients are not zero
loc<-coef_temp@i
##since the vector starts at position 1, so no need to use loc-1
marker<-rownames(train)[loc]
train.select<-train[rownames(train)%in%marker,]
#heatmap(t(train.select),labRow= train.classname,scale="none")



#http://sebastianraschka.com/Articles/heatmaps_in_r.html
if (!require("gplots")) {
  install.packages("gplots", dependencies = TRUE)
  library(gplots)
}
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", dependencies = TRUE)
  library(RColorBrewer)
}

# creates a own color palette from red to green
my_palette <- colorRampPalette(c("yellow", "red", "green"))(n = 299)
# (optional) defines the color breaks manually for a "skewed" color transition
col_breaks = c(seq(-1,0,length=100),  # for red
               seq(0.01,0.8,length=100),           # for yellow
               seq(0.81,1,length=100))             # for green
#matrix for heatmap

dm<-t(train.select)
rownames(dm)<-train.classname

distance = dist(dm, method = "manhattan")
cluster = hclust(distance, method = "ward.D2")
heatmap.2(dm,
          #cellnote = t(test.select),  # same data set for cell labels
          #main = "Clustering", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(4,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier
          breaks=col_breaks,    # enable color transition at specified limits
          dendrogram="row",     # only draw a row dendrogram
          Rowv=as.dendrogram(cluster),
          #Rowv = TRUE,
          scale="none",
          offsetRow=-0.2,
          offsetCol = -0.2,
          cexRow = 0.5,
          #cexCol = 0.5,
          key=FALSE
)




####Test the model
test.pos<-seq(1,length(class.list),1)[-train.pos]
test<-rm[,test.pos]
test<-data.preprocess(test)
pred<-predict(cvfit,newx = t(test),s="lambda.1se")

##Generate the predicted class
pred<-as.data.frame(pred)
predict.class<-rep(-1,nrow(pred))
predict.class[which(pred[,1]>0.5)]<-1
predict.class[which(abs(pred[,1])<0.5)]<- 0


##Generate the real class
real.class<-rep(0,nrow(pred))
test.pos<-seq(1,length(class.list),1)[-train.pos]
test_ALLT.pos<-which(test.pos%in% ALLT.pos)
test_ALLB.pos<-which(test.pos%in% ALLB.pos)
test_AML.pos<-which(test.pos%in% AML.pos)
real.class[test_ALLT.pos]<-(-1)
real.class[test_ALLB.pos]<-0
real.class[test_AML.pos]<-1

ALLT.accuracy<-1-length(which(real.class[test_ALLT.pos]!=predict.class[test_ALLT.pos]))/length(test_ALLT.pos)
ALLB.accuracy<-1-length(which(real.class[test_ALLB.pos]!=predict.class[test_ALLB.pos]))/length(test_ALLB.pos)
AML.accuracy<-1-length(which(real.class[test_AML.pos]!=predict.class[test_AML.pos]))/length(test_AML.pos)




accuracy<-1-length(which(real.class!=predict.class))/nrow(pred)
##check accuracy
print(c("ALL-T accuracy: ",ALLT.accuracy,"Test sample size: ",length(test_ALLT.pos)))
print(c("ALL-B accuracy: ",ALLB.accuracy,"Test sample size: ",length(test_ALLB.pos)))
print(c("AML accuracy: ",AML.accuracy,"Test sample size: ",length(test_AML.pos)))
print(c("Total accuracy: ",accuracy,"Total test sample size: ", length(test.pos)))



####Cluster the generated test dataset
test.pos<-seq(1,length(class.list),1)[-train.pos]
##pay attention! need to generate the relative positions in class.list, not just in the test.pos!
pos1<-ALLT.pos[which(ALLT.pos%in% test.pos )]
pos2<-ALLB.pos[which(ALLB.pos%in% test.pos )]
pos3<-AML.pos[which(AML.pos%in% test.pos )]

#class.list[pos1]

test.table<-getSample.pos(pos1,pos2,pos3,size = c(32,32,32),seed=940309)

test.pos<-test.table$pos
test.classlabel<-test.table$classlabel
test.classname<-class.list[test.pos]

##Generate the training dataset
test<-rm[,test.pos]

##Preprocess the training dataset
test<-data.preprocess(test)

##Draw heatmap for the selected test dataset by the markers generated from the training dataset
test.select<-test[rownames(test)%in%marker,]

#matrix for heatmap
dm<-t(test.select)
rownames(dm)<-test.classname

distance = dist(dm, method = "manhattan")
cluster = hclust(distance, method = "ward.D2")
heatmap.2(dm,
          #cellnote = t(test.select),  # same data set for cell labels
          #main = "Clustering", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(4,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier
          breaks=col_breaks,    # enable color transition at specified limits
          dendrogram="row",     # only draw a row dendrogram
          Rowv=as.dendrogram(cluster),
          #Rowv = TRUE,
          scale="none",
          offsetRow=-0.2,
          offsetCol = -0.2,
          cexRow = 0.5,
          #cexCol = 0.5,
          key=FALSE
)
