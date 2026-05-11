# load the required libraries
library("openxlsx")

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
mutations_chip = mutations[clinical_data$PATIENT_ID,]
mutations = NULL
load("analysis_data/MSK-ctDNA/final_data/mutations_2024.RData")
mutations = mutations_2024
mutations_2024 = NULL
valid_samples = sort(unique(intersect(clinical_data$PATIENT_ID,rownames(mutations))))
clinical_data = clinical_data[which(clinical_data$PATIENT_ID%in%valid_samples),]
rownames(clinical_data) = 1:nrow(clinical_data)
mutations_chip = mutations_chip[clinical_data$PATIENT_ID,]
colnames(mutations_chip) = paste0(colnames(mutations_chip),"_CHIP")
mutations_ctdna = mutations[clinical_data$PATIENT_ID,]
colnames(mutations_ctdna) = paste0(colnames(mutations_ctdna),"_ctDNA")
mutations = NULL
mutations = cbind(mutations_chip,mutations_ctdna)

# process the metastasis data
metastasis_data = clinical_data[,c("PATIENT_ID","AGE","GENDER","SAMPLE_TYPE","CANCER_TYPE_DETAILED")]
metastasis_data = metastasis_data[-which(is.na(metastasis_data),arr.ind=TRUE)[,"row"],]
metastasis_data$SAMPLE_TYPE[which(metastasis_data$SAMPLE_TYPE=="Local Recurrence")] = "Metastasis"
metastasis_data$GENDER[which(metastasis_data$GENDER=="Male")] = 0
metastasis_data$GENDER[which(metastasis_data$GENDER=="Female")] = 1
metastasis_data$GENDER = as.numeric(metastasis_data$GENDER)
cancer_types = "Lung Adenocarcinoma"
metastasis_data = metastasis_data[which(metastasis_data$CANCER_TYPE_DETAILED==cancer_types),]
rownames(metastasis_data) = 1:nrow(metastasis_data)
colnames(metastasis_data)[5] = "CANCER_TYPE"
mutations = mutations[metastasis_data$PATIENT_ID,]
mutations = mutations[,sort(unique(names(which(colSums(mutations)>=3))))]

# perform logistic regression with a binary target to associate mutations to metastasis
set.seed(12345)
model_estimates = list()
input_data = list()
for(i in sort(unique(metastasis_data$CANCER_TYPE))) {
    # process the data
    mets_data = metastasis_data[which(metastasis_data$CANCER_TYPE==i),]
    features_data = mutations[mets_data$PATIENT_ID,]
    features_data = features_data[,sort(unique(names(which((colSums(features_data)/nrow(features_data))>0.005)))),drop=FALSE]
    features_data = features_data[,sort(unique(names(which(colSums(features_data)>=3)))),drop=FALSE]
    mets_data$SAMPLE_TYPE[which(mets_data$SAMPLE_TYPE=="Primary")] = 0
    mets_data$SAMPLE_TYPE[which(mets_data$SAMPLE_TYPE=="Metastasis")] = 1
    mets_data$SAMPLE_TYPE = as.numeric(mets_data$SAMPLE_TYPE)
    analysis_data = data.frame(cbind(mets_data$SAMPLE_TYPE,mets_data$AGE,mets_data$GENDER,features_data))
    colnames(analysis_data)[c(1,2,3)] = c("Metastasis","AGE","GENDER")
    # perform multivariate logistic regression with a binary target
    model_logistic = glm(Metastasis ~ ., data = analysis_data, family = "binomial")
    model_logistic = summary(model_logistic)
    model_logistic = model_logistic$coefficients[,c("Estimate","Pr(>|z|)")]
    # save the results for the current cancer type
    model_estimate = data.frame(VARIABLE=gsub("\\.","-",rownames(model_logistic)),OR=exp(model_logistic[,"Estimate"]),PVALUE=model_logistic[,"Pr(>|z|)"])
    model_estimate = model_estimate[which(model_estimate$PVALUE<0.05),,drop=FALSE]
    if(length(which(rownames(model_estimate)=="(Intercept)"))>0) {
        model_estimate = model_estimate[-which(rownames(model_estimate)=="(Intercept)"),,drop=FALSE]
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
}

# save the results
save(input_data,file="results/input_data.RData")
save(model_estimates,file="results/model_estimates.RData")
write.xlsx(x=model_estimates,file="results/model_estimates.xlsx",rowNames=FALSE,colNames=TRUE)
