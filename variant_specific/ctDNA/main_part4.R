# load the required libraries
library("openxlsx")
library("survival")

# load the input data
load("analysis_data/MSK-SOLID/final_data/clinical_data_2021.RData")
load("analysis_data/MSK-SOLID/final_data/clinical_data_2024.RData")
clinical_data = rbind(clinical_data_2021,clinical_data_2024)
clinical_data = clinical_data[order(clinical_data$PATIENT_ID),]
clinical_data_2021 = NULL
clinical_data_2024 = NULL
load("analysis_data/MSK-ctDNA/final_data/mutations_2024.RData")
mutations = mutations_2024
mutations_2024 = NULL
valid_samples = sort(unique(intersect(clinical_data$PATIENT_ID,rownames(mutations))))
clinical_data = clinical_data[which(clinical_data$PATIENT_ID%in%valid_samples),]
rownames(clinical_data) = 1:nrow(clinical_data)
mutations = mutations[clinical_data$PATIENT_ID,]

# process the survival data
survival_data = clinical_data[,c("PATIENT_ID","AGE","GENDER","OS_MONTHS","OS_STATUS","CANCER_TYPE_DETAILED")]
colnames(survival_data) = c("PATIENT_ID","AGE","GENDER","OS_MONTHS","OS_STATUS","CANCER_TYPE_DETAILED")
survival_data = survival_data[which(!is.na(survival_data$AGE)),]
survival_data = survival_data[which(!is.na(survival_data$GENDER)),]
survival_data[,"OS_STATUS"][which(survival_data[,"OS_STATUS"]=="0:LIVING")] = 0
survival_data[,"OS_STATUS"][which(survival_data[,"OS_STATUS"]=="1:DECEASED")] = 1
survival_data[,"OS_STATUS"][which(as.numeric(survival_data[,"OS_MONTHS"])<1|as.numeric(survival_data[,"AGE"])>80)] = NA
survival_data[,"OS_MONTHS"][which(as.numeric(survival_data[,"OS_MONTHS"])<1|as.numeric(survival_data[,"AGE"])>80)] = NA
survival_data[,"OS_STATUS"][which(as.numeric(survival_data[,"OS_MONTHS"])>60)] = 0
survival_data[,"OS_MONTHS"][which(as.numeric(survival_data[,"OS_MONTHS"])>60)] = 60
survival_data$GENDER[which(survival_data$GENDER=="Male")] = 0
survival_data$GENDER[which(survival_data$GENDER=="Female")] = 1
survival_data$GENDER = as.numeric(survival_data$GENDER)
clinical_data = NULL
clinical_data = survival_data
rownames(clinical_data) = clinical_data$PATIENT_ID
invalid_samples = unique(which(is.na(clinical_data),arr.ind=TRUE)[,"row"])
if(length(invalid_samples)>0) {
    clinical_data = clinical_data[-invalid_samples,,drop=FALSE]
}
cancer_types = c("Breast Invasive Ductal Carcinoma","Lung Adenocarcinoma","Lung Squamous Cell Carcinoma","Pancreatic Adenocarcinoma","Prostate Adenocarcinoma")
clinical_data = clinical_data[which(clinical_data$CANCER_TYPE_DETAILED%in%cancer_types),]
rownames(clinical_data) = 1:nrow(clinical_data)
colnames(clinical_data)[6] = "CANCER_TYPE"
mutations = mutations[clinical_data$PATIENT_ID,]
mutations = mutations[,sort(unique(names(which(colSums(mutations)>=3))))]

# perform the survival association analysis
set.seed(12345)
model_estimates = list()
input_data = list()
for(i in sort(unique(clinical_data$CANCER_TYPE))) {
    # process the data
    survival_data = clinical_data[which(clinical_data$CANCER_TYPE==i),]
    rownames(survival_data) = survival_data$PATIENT_ID
    features_data = mutations[rownames(survival_data),]
    analysis_data = data.frame(cbind(survival_data[,c("OS_MONTHS","OS_STATUS","AGE","GENDER")],t(t(rowSums(features_data)))))
    colnames(analysis_data)[1:2] = c("Times","Status")
    colnames(analysis_data)[5] = "MUTATIONAL_BURDEN"
    # perform multivariate Cox regression
    analysis_cov = analysis_data[,colnames(analysis_data)[3:ncol(analysis_data)],drop=FALSE]
    string_test = paste0("analysis_cov$",colnames(analysis_cov),collapse="+")
    string_test = gsub("analysis_cov\\$","",string_test)
    time = as.numeric(analysis_data$Times)
    status = as.numeric(analysis_data$Status)
    string_test = paste0("Surv(time, status) ~ ",string_test,collapse="")
    string_test = as.formula(string_test)
    res_cox = coxph(formula = string_test, data = analysis_data)
    res_cox = summary(res_cox)
    res_cox = res_cox$coefficients[,c("exp(coef)","Pr(>|z|)")]
    # save the results for the current cancer type
    model_estimate = data.frame(VARIABLE=gsub("\\.","-",rownames(res_cox)),HR=res_cox[,"exp(coef)"],PVALUE=res_cox[,"Pr(>|z|)"])
    model_estimate = model_estimate[which(model_estimate$PVALUE<0.05),,drop=FALSE]
    if(nrow(model_estimate)>0) {
        input_data[[i]] = analysis_data
        rownames(model_estimate) = 1:nrow(model_estimate)
        model_estimates[[i]] = model_estimate
    }
}
selected_models = NULL
for(i in 1:length(model_estimates)) {
    selected_models = c(selected_models,length(grep("MUTATIONAL_BURDEN",model_estimates[[i]][,"VARIABLE"])))
}
model_estimates = model_estimates[which(selected_models!=0)]
names(model_estimates)[1] = "Breast Ductal Carcinoma"

# save the results
save(input_data,file="results/input_data_mb.RData")
save(model_estimates,file="results/model_estimates_mb.RData")
write.xlsx(x=model_estimates,file="results/model_estimates_mb.xlsx",rowNames=FALSE,colNames=TRUE)
