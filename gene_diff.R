library(xlsxjars)
library(rJava)
library(xlsx)
library(MASS)

#dataset downloaded from:
#https://archive.ics.uci.edu/ml/datasets/Mice+Protein+Expression


##first replace the NA by zero in microsoft excel(operation canceled)
raw_data<-read.csv("~/gene_diff/Data_Cortex_Nuclear.csv")


##delete the data of which the treatment is Memantine
data1<-raw_data[which(raw_data$Treatment != "Memantine"),]
##delete rows containing NA measurements
data1<-na.omit(data1)

#https://stackoverflow.com/questions/31139838/multiple-t-test-in-r
p_value_Genotype<-sapply(data1[,2:78], function(x) t.test(x ~ data1$Genotype)$p.value)

#https://www.rdocumentation.org/packages/stats/versions/3.4.0/topics/p.adjust
#also you can choose other method by the below command
#p.adjust.M<-p.adjust.methods
Bonferroni<-p.adjust(method="bonferroni",p_value_Genotype,n=77)

hub<-data.frame(Protein_or_modification,p_value_Genotype,Bonferroni)

hub_ranked<-hub[order(hub$Bonferroni),]
#marker<-rownames(hub)[which(hub$p_value_Genotype<0.001)]
marker<-rownames(hub)[which(hub$Bonferroni<0.001)]
##generate training dataset
included_cols<-c(marker,"Genotype")
#marker<-rownames(hub_ranked)[1:20]

##used to build a prediction model
data2<-data1[,which(colnames(data1)%in% included_cols)]



#https://stats.stackexchange.com/questions/61090/how-to-split-a-data-set-to-do-10-fold-cross-validation
#Randomly shuffle the data
temp_data<-data2[sample(nrow(data2)),]

#Create 10 equally size folds
folds <- cut(seq(1,nrow(temp_data)),breaks=10,labels=FALSE)

#Segement your data by fold using the which() function 
Indexes <- which(folds==1,arr.ind=TRUE)
test_data <- temp_data[Indexes, ]
train_data<-temp_data[-Indexes,]

##generate the lda model
train_lda<-lda(Genotype ~., data = train_data)
temp_df<-cbind(train_lda$scaling)

#Perform the rest 10 fold cross validation
for(i in 2:10){
  #Segement your data by fold using the which() function 
  Indexes <- which(folds==i,arr.ind=TRUE)
  test_data <- temp_data[Indexes, ]
  train_data <- temp_data[-Indexes, ]
  ##generate the lda model
  train_lda<-lda(Genotype ~., data = train_data)
  
  temp_df<-cbind(temp_df,train_lda$scaling)
}

temp_mtx<-as.matrix(temp_df)
coefficients<-rowMeans(temp_mtx)

?lda

train_lda$svd
pred<-predict(train_lda,test_data)$class
real<-test_data$Genotype
dim(train_data)
##check the wrong classfication rate
length(which(pred!=real))/nrow(test_data)


##generate the out_test_data as a control
out_test_data<-raw_data[which(raw_data$Treatment == "Memantine"),]
out_test_data<-na.omit(out_test_data)
out_test_data<-out_test_data[,which(colnames(out_test_data)%in% included_cols)]

pred<-predict(train_lda,out_test_data)$class
real<-out_test_data$Genotype
dim(train_data)
##check the wrong classification rate
length(which(pred!=real))/nrow(out_test_data)
