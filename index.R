# 사용자 정의 NRI 함수
calculate_nri <- function(event, old_probs, new_probs) {
  # 사건 발생 그룹과 사건 미발생 그룹으로 분류
  event_indices <- which(event == 1)
  nonevent_indices <- which(event == 0)
  
  # 사건 발생 그룹에서 상향 및 하향 조정 비율 계산
  event_upward_reclassification <- sum(new_probs[event_indices] > old_probs[event_indices]) / length(event_indices)
  event_downward_reclassification <- sum(new_probs[event_indices] < old_probs[event_indices]) / length(event_indices)
  
  # 사건 미발생 그룹에서 상향 및 하향 조정 비율 계산
  nonevent_upward_reclassification <- sum(new_probs[nonevent_indices] > old_probs[nonevent_indices]) / length(nonevent_indices)
  nonevent_downward_reclassification <- sum(new_probs[nonevent_indices] < old_probs[nonevent_indices]) / length(nonevent_indices)
  
  # NRI 계산
  nri <- (event_upward_reclassification - event_downward_reclassification) - (nonevent_upward_reclassification - nonevent_downward_reclassification)
  
  return(nri)
}

set.seed(42)

n_samples <- 100
old_model_probs <- runif(n_samples, min = 0, max = 1)
new_model_probs <- runif(n_samples, min = 0, max = 1)
true_labels <- rbinom(n_samples, size = 1, prob = 0.5)

# 실제 결과값을 이진형태로 변환 (1은 사건 발생, 0은 사건 미발생)
event <- as.numeric(true_labels == 1)

# NRI 계산
nri_result <- calculate_nri(event, old_model_probs, new_model_probs)

# 결과 출력
print(nri_result)


# 기존 모델과 새로운 모델의 예측 확률 및 실제 발생 여부
old_model_probs <- c(0.1, 0.4, 0.25, 0.7, 0.8)
new_model_probs <- c(0.2, 0.5, 0.35, 0.75, 0.9)
actual_outcomes <- c(0, 1, 0, 1, 1)

# 위험 범주의 경계값 정의
risk_thresholds <- c(0, 0.3, 0.6, 1)

# 기존 모델과 새로운 모델의 예측 결과를 위험 범주로 변환
old_model_categories <- findInterval(old_model_probs, risk_thresholds)
new_model_categories <- findInterval(new_model_probs, risk_thresholds)

# 사건 발생 그룹과 미발생 그룹으로 분리
event_group <- which(actual_outcomes == 1)
nonevent_group <- which(actual_outcomes == 0)

# 상향 조정 비율과 하향 조정 비율 계산
up_reclass_event <- sum(new_model_categories[event_group] > old_model_categories[event_group]) / length(event_group)
down_reclass_event <- sum(new_model_categories[event_group] < old_model_categories[event_group]) / length(event_group)
up_reclass_nonevent <- sum(new_model_categories[nonevent_group] > old_model_categories[nonevent_group]) / length(nonevent_group)
down_reclass_nonevent <- sum(new_model_categories[nonevent_group] < old_model_categories[nonevent_group]) / length(nonevent_group)

# Category NRI 계산
category_nri <- (up_reclass_event - down_reclass_event) - (up_reclass_nonevent - down_reclass_nonevent)
print(category_nri)

brier_score <- function(true_labels, predicted_probs) {
  mean((predicted_probs - true_labels)^2)
}


# 필요한 패키지 설치 및 불러오기
if (!requireNamespace("ROCR", quietly = TRUE)) {
  install.packages("ROCR")
}
library(ROCR)

# 모델이 예측한 확률 값을 사용하여 ROCR에서 제공하는 함수로 예측 객체 생성
pred <- prediction(predicted_probabilities, actual_labels)

# Precision-Recall Curve 계산
pr <- performance(pred, measure = "prec", x.measure = "rec")

# 그래프 그리기
plot(pr, main = "Precision-Recall Curve with Youden's J Threshold", col = "blue", lwd = 2)

# Youden's J statistic을 사용하여 구한 임계값에 대한 Precision과 Recall 표시
threshold <- y_threshold # 여기서 y_threshold는 Youden's J statistic을 사용하여 구한 임계값
binary_predictions <- ifelse(predicted_probabilities >= threshold, 1, 0)
precision <- posPredValue(binary_predictions, actual_labels)
recall <- sensitivity(binary_predictions, actual_labels)
points(recall, precision, col = "red", pch = 19)
text(recall, precision, labels = paste0("T = ", round(threshold, 2)), col = "red", pos = 4)

legend("bottomright", legend = c("Precision-Recall Curve", "Youden's J Threshold"), col = c("blue", "red"), pch = c(NA, 19), lty = c(1, NA), lwd = 2)

mod <- glm(as.factor(case) ~ scale(pred_inf) ,data=nri.data,family=binomial)
pred<- predict(mod,newdata=nri.data,type="response")

pred_obj <- ROCR::prediction(pred,nri.data$case)
roc_obj <- ROCR::performance(pred_obj,"sens","spec")
youden_index <- roc_obj@y.values[[1]] + roc_obj@x.values[[1]] - 1
optimal_threshold_index <- which.max(youden_index)
optimal_threshold <- roc_obj@alpha.values[[1]][optimal_threshold_index]

pr <- ifelse(pred >= optimal_threshold,1,0)
cm <- confusionMatrix(table(pr,actual_labels))
> cm
true_label <- nri.data$case
mod <- glm(as.factor(case) ~ scale(pred_inf),data=nri.data,family=binomial)
pred1 <- predict(mod,newdata=nri.data,type="response")
mod <- glm(as.factor(case) ~ scale(pred_inf)+ scale(metaPGS),data=nri.data,family=binomial)
pred2 <- predict(mod,newdata=nri.data,type="response")
roc_obj1 <- roc(true_label,pred1)
auc_obj1 <- auc(roc_obj1)
roc_obj2 <- roc(true_label,pred2)
auc_obj2 <- auc(roc_obj2)
plot(roc_obj1,main="ROC Curves Comparison",col="blue",lwd=2)
lines(roc_obj2,col="red",lwd=2)
abline(h=0:1,v=0:1,col="gray",lty=3)
legend("bottomright",legend=c("Model1","Model2"),col=c("blue","red"),lwd=2)
print(auc(roc_obj1))
Area under the curve: 0.6083
print(auc(roc_obj2))
Area under the curve: 0.6879


SMOTE <- function (X, target, K = 5, dup_size = 0) {
    ncD = ncol(X)
    n_target = as.numeric(table(target))
    classP = names(which.min(n_target))
    P_set = X[target == classP, ][sample(sum(target == classP)), ]
    N_set = X[target != classP, ]
    P_class = rep(classP, nrow(P_set))
    N_class = target[target != classP]
    sizeP = nrow(P_set)
    sizeN = nrow(N_set)
    knear = knearest(P_set, P_set, K)
    sum_dup = n_dup_max(sizeP + sizeN, sizeP, sizeN, dup_size)
    syn_dat = NULL
    for (i in 1:sizeP) {
        if (is.matrix(knear)) {
            pair_idx = knear[i, ceiling(runif(sum_dup) * K)]
        }
        else {
            pair_idx = rep(knear[i], sum_dup)
        }
        g = runif(sum_dup)
        P_i = matrix(unlist(P_set[i, ]), sum_dup, ncD, byrow = TRUE)
        Q_i = as.matrix(P_set[pair_idx, ])
        syn_i = P_i + g * (Q_i - P_i)
        syn_dat = rbind(syn_dat, syn_i)
    }
    P_set[, ncD + 1] = P_class
    colnames(P_set) = c(colnames(X), "class")
    N_set[, ncD + 1] = N_class
    colnames(N_set) = c(colnames(X), "class")
    rownames(syn_dat) = NULL
    syn_dat = data.frame(syn_dat)
    syn_dat[, ncD + 1] = rep(classP, nrow(syn_dat))
    colnames(syn_dat) = c(colnames(X), "class")
    NewD = rbind(P_set, syn_dat, N_set)
    rownames(NewD) = NULL
    D_result = list(data = NewD, syn_data = syn_dat, orig_N = N_set,
                    orig_P = P_set, K = K, K_all = NULL, dup_size = sum_dup,
                    outcast = NULL, eps = NULL, method = "SMOTE")
    class(D_result) = "gen_data"
    return(D_result)
}
function (X, target, K = 5, dup_size = 0) {
    ncD = ncol(X)
    n_target=table(target)
    #n_target = as.numeric(table(target))
    classP = as.integaer(names(which.min(n_target)))
    #P_set = X[target == classP, ][sample(sum(target == classP)), ]
    P_set <- X[target == classP, , drop = FALSE][sample(1:sum(target == classP)), , drop = FALSE]

    N_set = X[target != classP, ]
    N_set = X[target != classP, , drop=FALSE]
    P_class = rep(classP, nrow(P_set))
    N_class = target[target != classP]
    sizeP = nrow(P_set)
    sizeN = nrow(N_set)
    knear = knearest(P_set, P_set, K)
    sum_dup = n_dup_max(sizeP + sizeN, sizeP, sizeN, dup_size)
    syn_dat = NULL
    for (i in 1:sizeP) {
        if (is.matrix(knear)) {
            pair_idx = knear[i, ceiling(runif(sum_dup) * K)]
        }
        else {
            pair_idx = rep(knear[i], sum_dup)
        }
        g = runif(sum_dup)
        P_i = matrix(unlist(P_set[i, ]), sum_dup, ncD, byrow = TRUE)
        Q_i = as.matrix(P_set[pair_idx, ])
        syn_i = P_i + g * (Q_i - P_i)
        syn_dat = rbind(syn_dat, syn_i)
    }
    P_set[, ncD + 1] = P_class
    colnames(P_set) = c(colnames(X), "class")
    N_set[, ncD + 1] = N_class
    colnames(N_set) = c(colnames(X), "class")
    rownames(syn_dat) = NULL
    syn_dat = data.frame(syn_dat)
    syn_dat[, ncD + 1] = rep(classP, nrow(syn_dat))
    colnames(syn_dat) = c(colnames(X), "class")
    NewD = rbind(P_set, syn_dat, N_set)
    rownames(NewD) = NULL
    D_result = list(data = NewD, syn_data = syn_dat, orig_N = N_set,
                    orig_P = P_set, K = K, K_all = NULL, dup_size = sum_dup,
                    outcast = NULL, eps = NULL, method = "SMOTE")
    class(D_result) = "gen_data"
    return(D_result)
}
true_label <- nri.data$case
mod <- glm(as.factor(case) ~ scale(metaPGS),data=nri.data,family=binomial)
pred1 <- predict(mod,newdata=nri.data,type="response")
mod <- glm(as.factor(case) ~ scale(pred_inf),data=nri.data,family=binomial)
pred2 <- predict(mod,newdata=nri.data,type="response")
mod <- glm(as.factor(case) ~ scale(pred_inf)+ scale(metaPGS),data=nri.data,family=binomial)
pred3 <- predict(mod,newdata=nri.data,type="response")
roc_obj1 <- roc(true_label,pred1)
roc_obj2 <- roc(true_label,pred2)
roc_obj3 <- roc(true_label,pred3)

plot(roc_obj1,main="ROC Curves Comparison",col="blue",lwd=2)
lines(roc_obj2,col="red",lwd=2)
lines(roc_obj3,col="green",lwd=2)

print(auc(roc_obj1)
print(auc(roc_obj2))
print(auc(roc_obj3))
roc_ggplot1 <- ggroc(roc_obj1)
roc_ggplot2 <- ggroc(roc_obj2)
roc_ggplot3 <- ggroc(roc_obj3)

combined_roc_data <- rbind(roc_ggplot1$data,roc_ggplot2$data,roc_ggplot3$data)
  
roc_data1 <- data.frame(roc_obj1$sensitivities,roc_obj1$specificities,Model="Model1")
roc_data2 <- data.frame(roc_obj2$sensitivities,roc_obj2$specificities,Model="Model2")
roc_data3 <- data.frame(roc_obj3$sensitivities,roc_obj3$specificities,Model="Model3")
colnames(roc_data1) <- c("sensitivity", "specificity", "Model")
colnames(roc_data2) <- c("sensitivity", "specificity", "Model")
colnames(roc_data3) <- c("sensitivity", "specificity", "Model")

combined_roc_data <- rbind(roc_data1, roc_data2, roc_data3)
combined_roc_ggplot <- ggplot(data = combined_roc_data)

combined_roc_ggplot +
geom_line(aes(x = 1 - specificity, y = sensitivity, color = Model), size = 1) +
geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray", size = 1) +
labs(x = "1 - Specificity (False Positive Rate)", y = "Sensitivity (True Positive Rate)", title = "ROC Curves", color = "Legend")+
theme_minimal()

roc_data1$Model <- c("MetaPGS")
roc_data2$Model <- c("DiseasePGS")
roc_data3$Model <- c("Meta+DiseasePGS")

combined_roc_data <- rbind(roc_data1, roc_data2, roc_data3)
combined_roc_ggplot <- ggplot(data = combined_roc_data)

combined_roc_ggplot +
geom_line(aes(x = 1 - specificity, y = sensitivity, color = Model), size = 1) +
geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray", size = 1) +
labs(x = "1 - Specificity (False Positive Rate)", y = "Sensitivity (True Positive Rate)", title = "ROC Curves", color = "Legend")+
theme_minimal()





ifelse( i %in% names(pheno1), value <- pheno1[,c(1,which(i == names(pheno1)))],
ifelse( i %in% names(pheno2), value <- pheno2[,c(1,which(i == names(pheno2)))],
ifelse( i %in% names(pheno3), value <- pheno3[,c(1,which(i == names(pheno3)))],
ifelse( i %in% names(pheno4), value <- pheno4[,c(1,which(i == names(pheno4)))],
ifelse( i %in% names(pheno5), value <- pheno5[,c(1,which(i == names(pheno5)))],
value <- pheno6[,c(1,which(i == names(pheno6)))])))))


for (i in qt){
value <- read.csv(paste("../44_vali_pheno/",i,sep=""),sep="\t")
pgs <- read.csv(paste("../33_GPS_result/",i,sep=""),sep="\t")
tmp <- left_join(value,pgs,by="eid")
tmp2 <- na.omit(tmp)
i <- gsub("X","p",i)
i <- paste(i,"_i0",sep="")
ph <- which(i == names(tmp2))
print(cor.test(tmp2$pred_inf,tmp2[[ph]]))
}


gwas <- read.table(paste("./result.",i,".assoc.logistic",sep=""),header=T)
head(gwas)
str(gwas)
gwas <- na.omit(gwas)
head(gwas)
gwas <- left_join(gwas,root2,by="SNP")
head(gwas)
dim(gwas)
gwas2 <- na.omit(gwas)
dim(gwas2)
names(gwas2)
gwas2 <- gwas2[,c(1,2,3,13,14,15,8,12)]
gwas2$beta <-log(gwas2$OR)
names(gwas2)
gwas2 <- gwas2[,c(1,2,3,13,14,15,8,12)]
head(gwas2)
gwas2 <- gwas2[,c(1,2,3,5,4,6,7,8)]
head(gwas2)
names(gwas2) <- c("chr","rsid","pos","a0","a1","beta","beta_se","p")
head(gwas2)
pheno <- read.csv(paste("../input_plink/QC3.",i,sep=""),sep="\t")
head(pheno)
tab <- table(pheno$case)
tab
neff <- 4/((1/tab[[2]]) + (1/tab[[1]]))
neff
4/((1/15540) + (1/158948))
gwas2$n_eff <- neff
head(gwas2)
write.table(gwas2,paste("../input_GPS/",i,".csv",sep=""),sep=",",quote=F,row.names=F)
BT <- c('M17','M18','M20','M21','M23','M24','M25','M43','M46','M47','M48','M50','M51','M54','M65','M67','M70','M75','M76','M77','M79','M81','M86',
)
BT <- c('M20','M21','M23','M24','M25','M43','M46','M47','M48','M50','M51','M54','M65','M67','M70','M75','M76','M77','M79','M81','M86')
BT
for ( i in BT){
print(paste("disease is :",i,sep=""))
gwas <- read.table(paste("./result.",i,".assoc.logistic",sep=""),header=T)
print(str(gwas))
gwas <- na.omit(gwas)
gwas2 <- left_join(gwas,bim,by="SNP")
gwas2 <- na.omit(gwas2)
if( nrow(gwas) != nrow(gwas2)){
break
}
gwas2$beta <- log(gwas2$OR)
gwas2 <- gwas2[,c(1,2,3,13,14,15,8,12)]
gwas2 <- gwas2[,c(1,2,3,5,4,6,7,8)]
names(gwas2) <- c("chr","rsid","pos","a0","a1","beta","beta_se","p")
pheno <- read.csv(paste("../input_plink/QC3.",i,sep=""),sep="\t")
tab <- table(pheno$case)
neff <- 4/((1/tab[[2]]) + (1/tab[[1]]))
gwas2$n_eff <- neff
print(head(gwas2))
write.table(gwas2,paste("../input_GPS/",i,".csv",sep=""),sep=",",quote=F,row.names=F)
}
sumstats2 <- read.csv(paste("/BiO/Hyein/90Traits/BT/QT_BT/2nd_validation_GWAS/phase3_44_GWAS_BT_in_GWASset/input_GPS/",i,".csv",sep=""))
map <- obj.bigSNP$map[-(2:3)]
names(map) <- c("chr", "pos", "a0", "a1")
info_snp <- snp_match(sumstats2, map)
fam.order <- as.data.table(obj.bigSNP$fam)
setnames(fam.order, c("family.ID", "sample.ID"), c("FID", "IID"))
df_beta <- info_snp[,c("beta", "beta_se", "n_eff", "_NUM_ID_")]
ldsc <- snp_ldsc( ld, length(ld), chi2 = (df_beta$beta / df_beta$beta_se)^2, sample_size = df_beta$n_eff, blocks = NULL)
h2_est <- ldsc[["h2"]]
ind.test <- 1:nrow(G2)
beta_inf <- snp_ldpred2_inf(corr, df_beta, h2 = h2_est)
pred_inf <- big_prodVec( G2, beta_inf, ind.row = ind.test, ind.col = info_snp$`_NUM_ID_`)
names(fam.order)[1] <- c("FID")
fam.order <- cbind(fam.order,pred_inf)
vali_fam <- left_join(vali,fam.order,by="FID")
write.table(vali_fam,paste("/BiO/Hyein/90Traits/BT/QT_BT/2nd_validation_GWAS/phase3_44_GWAS_BT_in_GWASset/GPS/",i,sep=""),sep="\t",quote=F,row.names=F)


for ( i in adamBT){
sumstats2 <- read.csv(paste("/BiO/Hyein/90Traits/BT/QT_BT/2nd_validation_GWAS/phase3_44_GWAS_BT_in_GWASset/from_adam_input_GPS/",i,".csv",sep=""))
map <- obj.bigSNP$map[-(2:3)]
names(map) <- c("chr", "pos", "a0", "a1")
info_snp <- snp_match(sumstats2, map)
fam.order <- as.data.table(obj.bigSNP$fam)
setnames(fam.order, c("family.ID", "sample.ID"), c("FID", "IID"))
df_beta <- info_snp[,c("beta", "beta_se", "n_eff", "_NUM_ID_")]
ldsc <- snp_ldsc( ld, length(ld), chi2 = (df_beta$beta / df_beta$beta_se)^2, sample_size = df_beta$n_eff, blocks = NULL)
h2_est <- ldsc[["h2"]]
ind.test <- 1:nrow(G2)
beta_inf <- snp_ldpred2_inf(corr, df_beta, h2 = h2_est)
pred_inf <- big_prodVec( G2, beta_inf, ind.row = ind.test, ind.col = info_snp$`_NUM_ID_`)
names(fam.order)[1] <- c("FID")
fam.order <- cbind(fam.order,pred_inf)
vali_fam <- left_join(vali,fam.order,by="FID")
write.table(vali_fam,paste("/BiO/Hyein/90Traits/BT/QT_BT/2nd_validation_GWAS/phase3_44_GWAS_BT_in_GWASset/from_adam_GPS/",i,sep=""),sep="\t",quote=F,row.names=F)
}



for ( i in adamBT){
tryCatch({
sumstats2 <- read.csv(paste("/BiO/Hyein/90Traits/BT/QT_BT/2nd_validation_GWAS/phase3_44_GWAS_BT_in_GWASset/from_adam_input_GPS/",i,".csv",sep=""))
map <- obj.bigSNP$map[-(2:3)]
names(map) <- c("chr", "pos", "a0", "a1")
info_snp <- snp_match(sumstats2, map)
fam.order <- as.data.table(obj.bigSNP$fam)
setnames(fam.order, c("family.ID", "sample.ID"), c("FID", "IID"))
df_beta <- info_snp[,c("beta", "beta_se", "n_eff", "_NUM_ID_")]
ldsc <- snp_ldsc( ld, length(ld), chi2 = (df_beta$beta / df_beta$beta_se)^2, sample_size = df_beta$n_eff, blocks = NULL)
h2_est <- ldsc[["h2"]]
ind.test <- 1:nrow(G2)
beta_inf <- snp_ldpred2_inf(corr, df_beta, h2 = h2_est)
pred_inf <- big_prodVec( G2, beta_inf, ind.row = ind.test, ind.col = info_snp$`_NUM_ID_`)
names(fam.order)[1] <- c("FID")
fam.order <- cbind(fam.order,pred_inf)
vali_fam <- left_join(vali,fam.order,by="FID")
write.table(vali_fam,paste("/BiO/Hyein/90Traits/BT/QT_BT/2nd_validation_GWAS/phase3_44_GWAS_BT_in_GWASset/from_adam_GPS/",i,sep=""),sep="\t",quote=F,row.names=F)},error=function(e){cat("error disease is ",i,"\n",file = "/BiO/Hyein/90Traits/BT/QT_BT/2nd_validation_GWAS/phase3_44_GWAS_BT_in_GWASset/from_adam_GPS/GPS.log",append=T)})
}

for ( i in xBT){
tryCatch({
sumstats2 <- read.csv(paste("/BiO/Hyein/90Traits/BT/QT_BT/2nd_validation_GWAS/phase3_44_GWAS_BT_in_GWASset/from_x_input_GPS/",i,".csv",sep=""))
map <- obj.bigSNP$map[-(2:3)]
names(map) <- c("chr", "pos", "a0", "a1")
info_snp <- snp_match(sumstats2, map)
fam.order <- as.data.table(obj.bigSNP$fam)
setnames(fam.order, c("family.ID", "sample.ID"), c("FID", "IID"))
df_beta <- info_snp[,c("beta", "beta_se", "n_eff", "_NUM_ID_")]
ldsc <- snp_ldsc( ld, length(ld), chi2 = (df_beta$beta / df_beta$beta_se)^2, sample_size = df_beta$n_eff, blocks = NULL)
h2_est <- ldsc[["h2"]]
ind.test <- 1:nrow(G2)
beta_inf <- snp_ldpred2_inf(corr, df_beta, h2 = h2_est)
pred_inf <- big_prodVec( G2, beta_inf, ind.row = ind.test, ind.col = info_snp$`_NUM_ID_`)
names(fam.order)[1] <- c("FID")
fam.order <- cbind(fam.order,pred_inf)
vali_fam <- left_join(vali,fam.order,by="FID")
write.table(vali_fam,paste("/BiO/Hyein/90Traits/BT/QT_BT/2nd_validation_GWAS/phase3_44_GWAS_BT_in_GWASset/from_x_GPS/",i,sep=""),sep="\t",quote=F,row.names=F)},error=function(e){cat("error disease is ",i,"\n",file = "/BiO/Hyein/90Traits/BT/QT_BT/2nd_validation_GWAS/phase3_44_GWAS_BT_in_GWASset/from_x_GPS/GPS.log",append=T)})
}

for ( i in oldBT){
tryCatch({
sumstats2 <- read.csv(paste("/BiO/Hyein/90Traits/BT/QT_BT/2nd_validation_GWAS/phase3_44_GWAS_BT_in_GWASset/from_old_input_GPS/",i,".csv",sep=""))
map <- obj.bigSNP$map[-(2:3)]
names(map) <- c("chr", "pos", "a0", "a1")
info_snp <- snp_match(sumstats2, map)
fam.order <- as.data.table(obj.bigSNP$fam)
setnames(fam.order, c("family.ID", "sample.ID"), c("FID", "IID"))
df_beta <- info_snp[,c("beta", "beta_se", "n_eff", "_NUM_ID_")]
ldsc <- snp_ldsc( ld, length(ld), chi2 = (df_beta$beta / df_beta$beta_se)^2, sample_size = df_beta$n_eff, blocks = NULL)
h2_est <- ldsc[["h2"]]
ind.test <- 1:nrow(G2)
beta_inf <- snp_ldpred2_inf(corr, df_beta, h2 = h2_est)
pred_inf <- big_prodVec( G2, beta_inf, ind.row = ind.test, ind.col = info_snp$`_NUM_ID_`)
names(fam.order)[1] <- c("FID")
fam.order <- cbind(fam.order,pred_inf)
vali_fam <- left_join(vali,fam.order,by="FID")
write.table(vali_fam,paste("/BiO/Hyein/90Traits/BT/QT_BT/2nd_validation_GWAS/phase3_44_GWAS_BT_in_GWASset/from_old_GPS/",i,sep=""),sep="\t",quote=F,row.names=F)},error=function(e){cat("error disease is ",i,"\n",file = "/BiO/Hyein/90Traits/BT/QT_BT/2nd_validation_GWAS/phase3_44_GWAS_BT_in_GWASset/from_old_GPS/GPS.log",append=T)})
}

>>> for i in range(0,case.shape[0]):
...     dd = case.iloc[i]["p41270"]
...     dd = dd.split("|")
...     ind = []
...     for j in range(0,len(dd)):
...             if(m.match(dd[j])):
...                     ind.append(j)
...     dic[case.iloc[i]["eid"]] = ind


plot(1, 1, xlim = range(data2$onset), ylim = c(0, 0.2), type = "n",
xlab = "Age", ylab = "Cumulative incidence", main = "Cumulative incidence by group")


for ( g in unique(data2$scoresumtile4)){
data_group <- data2[data2$scoresumtile4 ==g,]
data_group <- data_group[order(data_group$onset),]
cuminc <- cumsum(data_group$case.x)/ nrow(data_group)
lines(data_group$onset,cuminc,col=group_colors[g],lty=g)
}
for ( i in 1:nrow(relDate)){
df <- relDate[i,"index"]
df <- c(strsplit(df,split=",")[[1]])
for ( j in df){
string = paste("p41280_a",j,sep="")
tmp = relDate[i,string]
tmp = c(strsplit(tmp,split="-")[[1]])[1]
tmp <- as.integer(tmp)
min <- tmp
if(min > tmp){
min <- tmp
}
}
relDate[i,"inci"] <- min
}
