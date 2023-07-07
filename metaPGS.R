library(dplyr)
library(glmnet)
#########################setting###############################




#get association
assoc <- read.csv("asso.txt",sep="\t")
disease <- assoc$ICD10
dd <- unique(disease)
disease <- dd


#get risk factors with a single disease for unrelatedsample
tmp <- c()
for ( i in disease){
df <- assoc[assoc$ICD10 == i,]
tmp <- df$Pheno
assign(paste(i,"rf",sep=""),tmp)
}

#set unrelated sample root 
unrelated_root <- read.csv("/BiO/Hyein/90Traits/BT/QT_BT/2nd_validation_GWAS/11_rawQC/22_vali_root.csv")
df <- unrelated_root[,c(1,2,18,3,4,6:16)]
unrelated_root <- df


#calling all selected risk factors PGS for unrelatedsample
riskfactor <- assoc$Pheno
riskfactor <- unique(riskfactor)
#115 risk factor in unrelatedsample
for ( i in riskfactor){
assign(paste("PGS",i,sep=""),read.csv(paste("/BiO/Hyein/90Traits/BT/QT_BT/2nd_validation_GWAS/44_GPS_outcome/X",i,sep="")))
}
for ( i in riskfactor){
tmp <- get(paste("PGS",i,sep=""))
tmp <- tmp[,c(1,7)]
names(tmp)[2] <- c(paste("PGS",i,sep=""))
assign(paste("PGS",i,"v2",sep=""),tmp)
}
for( i in riskfactor){
tmp <- get(paste("PGS",i,"v2",sep=""))
sc <- scale(tmp[[2]])
nm <- c(paste("st_PGS",i,sep=""))
tmp[nm] <- sc
assign(paste("PGS",i,"v2",sep=""),tmp)
}
for ( i in riskfactor){
df <- get(paste("PGS",i,"v2",sep=""))
df <- df[,c(1,3)]
assign(paste("PGS",i,"v2",sep=""),df)
}

#get Binary Trait(BT;Disease) PGS in unrelatedsample(170,000)
for ( i in disease){
 assign(paste("PGS",i,sep=""),read.csv(paste("/BiO/Hyein/90Traits/BT/QT_BT/2nd_validation_GWAS/phase3_44_GWAS_BT_in_GWASset/outcome_GPS_GPSset/",i,sep="")))
}
for ( i in disease){
 tmp <- get(paste("PGS",i,sep=""))
 tmp <- tmp[,c("FID","pred_inf")]
 names(tmp)[2] <- c(paste("PGS",i,sep=""))
 sc  <- scale(tmp[[2]])
 nm <- c(paste("st_PGS",i,sep=""))
 tmp[nm] <- sc
tmp <-  tmp[,c(1,3)]
 assign(paste("PGS",i,"v2",sep=""),tmp)
 }

#paste a single disease case/control
for( i in disease){
if( i %in% commonDisease){
assign(paste(i,"case",sep=""),read.csv(paste("/BiO/Hyein/90Traits/BT/QT_BT/2nd_validation_GWAS/66_ICD10/22_extract_ICD10_sample_QC/",i,"/case.txt",sep=""),sep="\t"))
df <- get(paste(i,"case",sep=""))
df$case <- 1
df <- df[,c("eid","case")]
assign(paste(i,"case",sep=""),df)
}else{
assign(paste(i,"case",sep=""),read.csv(paste("/BiO/Hyein/90Traits/BT/QT_BT/2nd_validation_GWAS/66_ICD10/22_extract_ICD10_sample_QC/",i,sep=""),sep="\t",header=F))
df <- get(paste(i,"case",sep=""))
names(df)[1] <- c("eid")
df$case <- 1
df <- df[,c("eid","case")]
assign(paste(i,"case",sep=""),df)
}
}


#################Elasticnet############################
for ( i in disease){
unrmodel <- unrelated_root

for ( j in get(paste(i,"rf",sep=""))){
df <- get(paste("PGS",j,"v2",sep=""))
unrmodel <- left_join(unrmodel,df,by="FID")
}
df <- get(paste("PGS",i,"v2",sep=""))
unrmodel <- unrmodel(unrmodel,df,by="FID")
df <- get(paste(i,"case",sep=""))
unrmodel <- unrmodel(unrmodel,df,by="eid")
unrmodel[is.na(unrmodel$case),]$case <- 0

#perform Elasticnet logistic regression
train <- unrmodel[,c(-1,-2,-3)]
#train <- train[,c(14:length(names(train)))]
train.x <- scale(train[,-length(names(train))])
train.y <- train[,length(names(train))]
Elasticnet_model <- glmnet(train.x, train.y, alpha=0.5)
set.seed(100)
cv.Elasticnet <- cv.glmnet(train.x,train.y,alpha=0.5)
(best_lambda <- cv.Elasticnet$lambda.min)
(Elasticnet_coef <- coef(Elasticnet_model,s=best_lambda))

#select risk facotrs and opt lambda from elasticnet logisticregression 
ind <- Elasticnet_coef@i
rf <- Elasticnet_coef@Dimnames[[1]]
tmp <- c()
for ( j in 1:length(ind)){
tmp <- c(tmp,rf[ind[j]+1])
}
optrf <- tmp

#save elasticnet logistic regression result
assign(paste(i,"ElasticnetResult",sep=""),data.frame(risk.factor=optrf,opt_lambda=Elasticnet_coef@x))
df <- get(paste(i,"ElasticnetResult",sep=""))
write.table(df,paste("/BiO/Hyein/90Traits/BT/QT_BT/2nd_validation_GWAS/phase3_33_Elasticnet/QTBT/",i,"Result.txt",sep=""),sep="\t",quote=F,row.names=F)
}
###########Consctruction of QTmetPRS,QTBTmetaPRS####################
#QTmetaPRS
for ( d in sigdisease){
assign(paste("X",d,"Result.v1",sep=""),read.csv(paste("./QT/",d,"Result.txt",sep=""),sep="\t"))
df <- get(paste("X",d,"Result.v1",sep=""))
resultrf <- df$risk.factor
ind <- grep("st_PGS",resultrf)
resultrf <- resultrf[ind]
resultrf <- gsub("st_PGS","",resultrf)

reweight <- info_snp[,c(5,3)]
for ( i in resultrf){
df <- get(paste("X",i,".beta.inf",sep=""))
df <- df[,c(5,12)]
weight <- get(paste("X",d,"Result.v1",sep=""))
str <- c(paste("st_PGS",i,sep=""))
w <- weight[weight$risk.factor == str,]$opt_lambda
r <- SD[SD$riskfactor ==i,]$PGSsetSD
df$reweight <- df$beta_inf * (w/r)
df <- df[,c(1,3)]
names(df)[2] <- c(paste("X",i,".reweight.beta.inf",sep=""))
reweight <- left_join(reweight,df,by="rsid")
}

for ( i in 3:length(names(reweight))){
if(i == 3){
reweight$sum <- reweight[[i]]
}
else{
reweight$sum <- reweight$sum + reweight[[i]]
}
}

input <- reweight[,c("rsid","a0","sum")]
names(input) <- c("rsID","effect_allele","effect_weight")

if(!dir.exists(paste("./perSNPscore/QTSNPweight/",d,"/",sep=""))){
dir.create(paste("./perSNPscore/QTSNPweight/",d,"/",sep=""))
}
write.table(input,paste("./perSNPscore/QTSNPweight/",d,"/",d,".v1.txt",sep=""),sep="\t",quote=F,row.names=F)


}



#QTBTmetaPGS
for ( d in sigdisease){


assign(paste("X",d,"Result.v2",sep=""),read.csv(paste("./QTBT/",d,"Result.v2.txt",sep=""),sep="\t"))
df <- get(paste("X",d,"Result.v2",sep=""))
resultrf <- df$risk.factor
ind <- grep("st_PGS",resultrf)
resultrf <- resultrf[ind]
resultrf <- gsub("st_PGS","",resultrf)

reweight <- info_snp[,c(5,3)]
for ( i in resultrf){
df <- get(paste("X",i,".beta.inf",sep=""))
df <- df[,c(5,12)]
weight <- get(paste("X",d,"Result.v2",sep=""))
str <- c(paste("st_PGS",i,sep=""))
w <- weight[weight$risk.factor == str,]$opt_lambda
r <- SD[SD$riskfactor ==i,]$PGSsetSD
df$reweight <- df$beta_inf * (w/r)
df <- df[,c(1,3)]
names(df)[2] <- c(paste("X",i,".reweight.beta.inf",sep=""))
reweight <- left_join(reweight,df,by="rsid")
}

for ( i in 3:length(names(reweight))){
if(i == 3){
reweight$sum <- reweight[[i]]
}
else{
reweight$sum <- reweight$sum + reweight[[i]]
}
}

input <- reweight[,c("rsid","a0","sum")]
names(input) <- c("rsID","effect_allele","effect_weight")
if(!dir.exists(paste("./perSNPscore/QTBTSNPweight/",d,"/",sep=""))){
dir.create(paste("./perSNPscore/QTBTSNPweight/",d,"/",sep=""))
}
write.table(input,paste("./perSNPscore/QTBTSNPweight/",d,"/",d,".v2.txt",sep=""),sep="\t",quote=F,row.names=F)

}


#########QTmetaPGS,QTBTmetaPGS,BTPGS merge############################
for ( d in sigdisease){
v1 <- read.table(paste("/BiO/Hyein/90Traits/BT/QT_BT/2nd_validation_GWAS/phase3_33_Elasticnet/perSNPscore/QTSNPweight/",d,"/",d,".v1.profile",sep=""),header=T)
v2 <- read.table(paste("/BiO/Hyein/90Traits/BT/QT_BT/2nd_validation_GWAS/phase3_33_Elasticnet/perSNPscore/QTBTSNPweight/",d,"/",d,".v2.profile",sep=""),header=T)
bt <- read.csv(paste("/BiO/Hyein/90Traits/BT/QT_BT/2nd_validation_GWAS/phase3_44_GWAS_BT_in_GWASset/outcome_GPS_Repset/",d,sep=""))

v1 <- v1[,c("FID","SCORESUM")]
v2 <- v2[,c("FID","SCORESUM")]
bt <- bt[,c("FID","pred_inf")]

names(v1)[2] <- c("QTmetaPGS")
names(v2)[2] <- c("QTBTmetaPGS")
names(bt)[2] <- c("BTPGS")

merge <- left_join(v1,v2,by="FID")
merge <- left_join(merge,bt,by="FID")

merge <- left_join(merge,id,by="FID")
merge <- left_join(merge,cov,by="eid")
df <- get(paste(d,"case",sep=""))
merge <- left_join(merge,df,by="eid")
merge[is.na(merge$case),]$case <- 0
if(!dir.exists(paste("/BiO/Hyein/90Traits/BT/QT_BT/2nd_validation_GWAS/phase3_33_Elasticnet/perSNPscore/diseaseResult/",d,sep=""))){
 dir.create(paste("/BiO/Hyein/90Traits/BT/QT_BT/2nd_validation_GWAS/phase3_33_Elasticnet/perSNPscore/diseaseResult/",d,sep=""))
}
write.table(merge,paste("/BiO/Hyein/90Traits/BT/QT_BT/2nd_validation_GWAS/phase3_33_Elasticnet/perSNPscore/diseaseResult/",d,"/",d,".merge.txt",sep=""),sep="\t",quote=F,row.names=F)
 }
 










