library(dplyr)
library(glmnet)

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

for ( i in disease){
#i <- disease[i]
#paste a single disease case 
unrmodel <- unrelated_root

for ( j in get(paste(i,"rf",sep=""))){
df <- get(paste("PGS",j,"v2",sep=""))
unrmodel <- left_join(unrmodel,df,by="FID")
}
#paste a single disease case/control
if( i %in% cd){
assign(paste(i,"case",sep=""),read.csv(paste("/BiO/Hyein/90Traits/BT/QT_BT/2nd_validation_GWAS/66_ICD10/22_extract_ICD10_sample_QC/",i,"/case.txt",sep=""),sep="\t"))
}else{
assign(paste(i,"case",sep=""),read.csv(paste("/BiO/Hyein/90Traits/BT/QT_BT/2nd_validation_GWAS/66_ICD10/22_extract_ICD10_sample_QC/",i,sep=""),sep="\t",header=F))
}
df <- get(paste(i,"case",sep=""))
if( length(colnames(df)) > 1){
df <- df[,c(1)]
df <- data.frame(df)
names(df)[1] <- c("eid")
assign(paste(i,"case",sep=""),df)
}
names(df)[1] <- c("eid")
df["case"] <- 1
unrmodel <- left_join(unrmodel,df,by="eid")
assign(paste(i,"case",sep=""),df)
unrmodel[is.na(unrmodel$case),]$case <- 0

#perform Elasticnet logistic regression
train <- unrmodel[,c(-1,-2,-3)]
train <- train[,c(14:length(names(train)))]
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
write.table(df,paste("./",i,"Result.txt",sep=""),sep="\t",quote=F,row.names=F)

#For calculate metaPGS, call 50,000samples 
rel_root <- read.csv("/BiO/Hyein/90Traits/BT/QT_BT/2nd_validation_GWAS/phase2_11_related_sample/00_keep/11_related_sample_final.csv")
releval <- rel_root[,c(1,8)]
cov <- read.csv("/BiO/Hyein/90Traits/BT/QT_BT/2nd_validation_GWAS/phase2_11_related_sample/00_keep/ALL_cov.txt",sep="\t")
releval <- left_join(releval,cov,by="eid")
releval <- releval[,c(1:4,6:15,18)]
names(releval)[length(names(releval))] <- c("array")
df <- get(paste(i,"case",sep=""))
releval <- left_join(releval,df,by="eid")
releval[is.na(releval$case),]$case <- 0

#paste elasticnet optlambda risk factor PGS in 50,000 samples to the single disease
for ( j in get(paste(i,"rf",sep=""))){
df <- get(paste("relPGS",j,"v2",sep=""))
releval <- left_join(releval,df,by="eid")
}
df <- get(paste(i,"ElasticnetResult",sep=""))
#dd <- df[df$risk.factor != "(Intercept)" & df$risk.factor != "sex" & df$risk.factor != "age" & 
 #        df$risk.factor != "array" & df$risk.factor != "pc1" & df$risk.factor != "pc2" & 
  #       df$risk.factor != "pc3" & df$risk.factor != "pc4" & df$risk.factor != "pc5" 
  #       & df$risk.factor != "pc6" & df$risk.factor != "pc7" & df$risk.factor != "pc8" 
  #       & df$risk.factor != "pc9" & df$risk.factor != "pc10",]
Incp <- 0
Incp <- df[df$risk.factor == "(Intercept)","opt_lambda"]
dd <- df[!df$risk.factor == "(Intercept)",]
optrf <- dd$risk.factor
optlambda <- dd$opt_lambda
tmp <- gsub("st_PGS","st_relPGS",optrf)

#calculate metaPGS excluding covariated
relevalopt <- releval[,c(tmp,"case","eid","FID")]
coln <- names(relevalopt)
dd$risk.factor <- tmp
for ( j in tmp){
p <- which(j == coln)
Incp <- Incp + (( dd[dd$risk.factor == j,]$opt_lambda) * relevalopt[[p]])
}
relevalopt["metaPGS"] <- Incp

mod <- glm(as.factor(case) ~ scale(metaPGS),data=relevalopt,family=binomial)
allbeta <- summary(glm(as.factor(case) ~ scale(metaPGS),data=relevalopt,family=binomial))$coef["scale(metaPGS)",1]
allse <- summary(glm(as.factor(case) ~ scale(metaPGS),data=relevalopt,family=binomial))$coef["scale(metaPGS)",2]
allPvalue <- summary(glm(as.factor(case) ~ scale(metaPGS),data=relevalopt,family=binomial))$coef["scale(metaPGS)",4]
relevalopt <- relevalopt %>% mutate(tile100 = ntile(metaPGS,100))
write.table(relevalopt,paste("./",i,"metaPGS.txt",sep=""),sep="\t",quote=F,row.names=F)
nkr <- NagelkerkeR2(mod)$R2
g3 <- relevalopt[relevalopt$tile100 > 97,]
g3["marker"] <- 1
g4060 <- relevalopt[relevalopt$tile100 > 40,]
g4060 <- g4060[g4060$tile100 < 61,]
g4060["marker"] <- 0
merge <- rbind(g3,g4060)
allbeta <- summary(glm(as.factor(case) ~ scale(metaPGS),data=relevalopt,family=binomial))$coef["scale(metaPGS)",1]
allse <- summary(glm(as.factor(case) ~ scale(metaPGS),data=relevalopt,family=binomial))$coef["scale(metaPGS)",2]
allPvalue <- summary(glm(as.factor(case) ~ scale(metaPGS),data=relevalopt,family=binomial))$coef["scale(metaPGS)",4]

HRbeta <- summary(glm(as.factor(case) ~ as.factor(marker),data=merge,family=binomial))$coef["as.factor(marker)1",1]
HRse <- summary(glm(as.factor(case) ~ as.factor(marker),data=merge,family=binomial))$coef["as.factor(marker)1",2]
HRPvalue <- summary(glm(as.factor(case) ~ as.factor(marker),data=merge,family=binomial))$coef["as.factor(marker)1",4]

cat(i,"metaPGSallOR:",exp(allbeta)," ",i,"metaPGSallbeta:",allbeta," ",i,"metaPGSallse:",allse," ",i,"metaPGSallPvalue:",allPvalue,"\n",file="./output.txt",append=T)
cat(i,"metaPGSallNagelkerkeR2:",nkr,"\n",file="./output.txt",append=T)
cat(i,"metaPGSHROR:",exp(HRbeta)," ",i,"metaPGSHRbeta:",HRbeta," ",i,"metaPGSHRse:",HRse," ",i,"metaPGSHRPvalue:",HRPvalue,"\n",file="./output.txt",append=T)
}
