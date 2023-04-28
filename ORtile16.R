library(dplyr)

$diseasemetaPGS <- read.csv("/BiO/Hyein/90Traits/BT/QT_BT/2nd_validation_GWAS/phase3_33_Elasticnet/$diseasemetaPGS.txt",sep="\t")
$diseasePGS <- read.csv("/BiO/Hyein/",sep="\t")

diseasemetaPGS <- diseasemetaPGS[,c("eid","FID","case","metaPGS")]
diseasePGS <- diseasePGS[,c("eid","FID","pred_inf")]
names(diseasePGS)[3] <- c("diseasePGS")

merge <- left_join(diseasemetaPGS,diseasePGS,by="eid")

merge <- merge %>% mutate(diseasemetaPGStile4 = ntile(metaPGS,4))
merge <- merge %>% mutate(diseasePGStile4 = ntile(diseasePGS,4))

ref <- merge[merge$diseasePGStile4 == 1 & merge$diseasemetaPGStile4 ==1,]

for ( i in 1:4){
if ( i ==1){
for ( j in 2:4){
tmp <- merege[merge$diseasePGStile4 == j & merge$diseasemetaPGStile4 == i,]
df <- rbind(ref,tmp)
df$marker <- ifelse(df$diseasePGStile4 == 1 & df$diseasemetaPGStile4 == 1,0,1)
beta <- summary(glm(as.factor(case) ~ as.factor(marker),data=df,family=binomial))$coef["as.factor(marker)1",1]
Pvalue <- summary(glm(as.factor(case) ~ as.factor(marker),data=df,family=binomial))$coef["as.factor(marker)1",4]
print(paste("Beta: ",beta," OR: ",exp(beta)," Pvalue: ",Pvalue,sep=""))
}
else{
for ( j in 1:4){
tmp <- merge[merge$diseasePGStile4 == j & merge$diseasemetaPGStile4 == i,]
df <- rbind(ref,tmp)
df$marker <- ifelse(df$diseasePGStile4 == 1 & df$diseasemetaPGStile4 ==1,0,1)
beta <- summary(glm(as.factor(case) ~ as.factor(marker),data=df,family=binomial))$coef["as.factor(marker)1",1] 
Pvalue <- summary(glm(as.factor(case) ~ as.factor(marker),data=df,family=binomial))$coef["as.factor(marker)1",4]
print(paste("Beta: ",beta," OR: ",exp(beta)," Pvalue: ",Pvalue,sep=""))
}
}
}
