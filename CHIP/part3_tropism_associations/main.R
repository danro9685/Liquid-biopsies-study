# load the required libraries
library("openxlsx")
library("VGAM")

# load the input data
load("analysis_data/MSK-SOLID/final_data/clinical_data_2021.RData")
load("analysis_data/MSK-SOLID/final_data/clinical_data_2024.RData")
clinical_data = rbind(clinical_data_2021,clinical_data_2024)
clinical_data = clinical_data[order(clinical_data$PATIENT_ID),]
clinical_data_2021 = NULL
clinical_data_2024 = NULL
load("analysis_data/MSK-CHIP/final_data/mutations_2023.RData")
mutations = mutations_2023
mutations_2023 = NULL
valid_samples = sort(unique(intersect(clinical_data$PATIENT_ID,rownames(mutations))))
clinical_data = clinical_data[which(clinical_data$PATIENT_ID%in%valid_samples),]
rownames(clinical_data) = 1:nrow(clinical_data)
mutations = mutations[clinical_data$PATIENT_ID,]

# process the metastasis data
metastasis_data = clinical_data[,c("PATIENT_ID","AGE","GENDER","SAMPLE_TYPE","METASTATIC_SITE","CANCER_TYPE_DETAILED")]
metastasis_data$SAMPLE_TYPE[which(metastasis_data$SAMPLE_TYPE=="Local Recurrence")] = "Metastasis"
metastasis_data$METASTATIC_SITE[which(metastasis_data$SAMPLE_TYPE=="Primary")] = "Primary"
metastasis_data$SAMPLE_TYPE = NULL
metastasis_data$GENDER[which(metastasis_data$GENDER=="Male")] = 0
metastasis_data$GENDER[which(metastasis_data$GENDER=="Female")] = 1
metastasis_data$GENDER = as.numeric(metastasis_data$GENDER)
metastasis_data = metastasis_data[-which(is.na(metastasis_data),arr.ind=TRUE)[,"row"],]
cancer_types = c("Bladder Urothelial Carcinoma","Breast Invasive Ductal Carcinoma","Breast Invasive Lobular Carcinoma","Colon Adenocarcinoma","Cutaneous Melanoma","Esophageal Adenocarcinoma","Gastrointestinal Stromal Tumor","Hepatocellular Carcinoma","High-Grade Serous Ovarian Cancer","Intrahepatic Cholangiocarcinoma","Lung Adenocarcinoma","Lung Squamous Cell Carcinoma","Pancreatic Adenocarcinoma","Pancreatic Neuroendocrine Tumor","Papillary Thyroid Cancer","Prostate Adenocarcinoma","Rectal Adenocarcinoma","Renal Clear Cell Carcinoma","Small Cell Lung Cancer","Stomach Adenocarcinoma","Upper Tract Urothelial Carcinoma","Uterine Carcinosarcoma/Uterine Malignant Mixed Mullerian Tumor","Uterine Endometrioid Carcinoma","Uterine Serous Carcinoma/Uterine Papillary Serous Carcinoma")
metastasis_data = metastasis_data[which(metastasis_data$CANCER_TYPE_DETAILED%in%cancer_types),]
rownames(metastasis_data) = 1:nrow(metastasis_data)
colnames(metastasis_data)[5] = "CANCER_TYPE"
metastasis_data$CANCER_TYPE[which(metastasis_data$CANCER_TYPE=="Uterine Carcinosarcoma/Uterine Malignant Mixed Mullerian Tumor")] = "Uterine Carcinosarcoma"
metastasis_data$CANCER_TYPE[which(metastasis_data$CANCER_TYPE=="Uterine Serous Carcinoma/Uterine Papillary Serous Carcinoma")] = "Uterine Serous Carcinoma"
mutations = mutations[metastasis_data$PATIENT_ID,]
mutations = mutations[,sort(unique(names(which(colSums(mutations)>=3))))]

# perform multinomial logistic regression with a chategorical target to associate mutations to metastatic sites
set.seed(12345)
model_estimates = list()
input_data = list()
for(i in sort(unique(metastasis_data$CANCER_TYPE))[-20]) {
    # process the data
    mets_data = metastasis_data[which(metastasis_data$CANCER_TYPE==i),]
    features_data = mutations[mets_data$PATIENT_ID,]
    features_data = features_data[,sort(unique(names(which((colSums(features_data)/nrow(features_data))>0.005)))),drop=FALSE]
    features_data = features_data[,sort(unique(names(which(colSums(features_data)>=3)))),drop=FALSE]
    mets_data$METASTATIC_SITE[which(mets_data$METASTATIC_SITE%in%names(which(table(mets_data$METASTATIC_SITE)<5)))] = "Other"
    analysis_data = data.frame(cbind(mets_data$METASTATIC_SITE,mets_data$AGE,mets_data$GENDER,features_data))
    colnames(analysis_data)[c(1,2,3)] = c("Metastasis","AGE","GENDER")
    analysis_data$Metastasis = factor(analysis_data$Metastasis)
    analysis_data$AGE = as.numeric(analysis_data$AGE)
    analysis_data$GENDER = as.numeric(analysis_data$GENDER)
    for(j in 4:ncol(analysis_data)) {
        analysis_data[,j] = as.numeric(analysis_data[,j])
    }
    constant_cols = sapply(analysis_data, function(x) length(unique(x)) == 1)
    analysis_data = analysis_data[,!constant_cols]
    # perform multinomial logistic regression with a chategorical target
    model_logistic = vglm(Metastasis ~ ., data = analysis_data, family = "multinomial", refLevel = "Primary")
    model_logistic = summary(model_logistic)
    model_logistic = coef(model_logistic)
    model_logistic = model_logistic[,c("Estimate","Pr(>|z|)")]
    # save the results for the current cancer type
    sites_no_ref = setdiff(levels(analysis_data$Metastasis),"Primary")
    covariates = gsub("\\.","-",rownames(model_logistic))
    VARIABLE = NULL
    SITE = NULL
    for(j in covariates) {
        VARIABLE = c(VARIABLE,strsplit(j,split=":")[[1]][[1]])
        SITE = c(SITE,sites_no_ref[as.numeric(strsplit(j,split=":")[[1]][[2]])])
    }
    model_estimate = data.frame(VARIABLE=VARIABLE,SITE=SITE,OR=exp(model_logistic[,"Estimate"]),PVALUE=model_logistic[,"Pr(>|z|)"])
    model_estimate = model_estimate[which(model_estimate$PVALUE<0.05),,drop=FALSE]
    if(length(which(model_estimate$VARIABLE=="(Intercept)"))>0) {
        model_estimate = model_estimate[-which(model_estimate$VARIABLE=="(Intercept)"),,drop=FALSE]
    }
    if(nrow(model_estimate)>0) {
        input_data[[i]] = analysis_data
        rownames(model_estimate) = 1:nrow(model_estimate)
        NUM_MUT = rep(NA,nrow(model_estimate))
        for(j in 1:nrow(model_estimate)) {
            if(model_estimate$VARIABLE[j]!="AGE"&&model_estimate$VARIABLE[j]!="GENDER") {
                NUM_MUT[j] = sum(features_data[,model_estimate$VARIABLE[j]])
            }
        }
        FREQ_MUT = (NUM_MUT/nrow(features_data))
        model_estimate$NUM_MUT = NUM_MUT
        model_estimate$FREQ_MUT = FREQ_MUT
        model_estimates[[i]] = model_estimate
    }
    cat(i,"\n")
}
names(model_estimates)[2] = "Breast Ductal Carcinoma"
names(model_estimates)[3] = "Breast Lobular Carcinoma"
names(model_estimates)[7] = "High-Grade Ovarian Cancer"

# save the results
save(input_data,file="results/input_data.RData")
save(model_estimates,file="results/model_estimates.RData")
write.xlsx(x=model_estimates,file="results/model_estimates.xlsx",rowNames=FALSE,colNames=TRUE)
