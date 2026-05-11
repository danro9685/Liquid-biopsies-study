# load the required libraries
library("openxlsx")
library("survival")

# load and process the input data
load("results/clinical_data.RData")
load("results/mutations_data.RData")
clinical_data = clinical_data[which(clinical_data$PATIENT_ID%in%intersect(clinical_data$PATIENT_ID,mutations_data$PATIENT_ID)),]
rownames(clinical_data) = 1:nrow(clinical_data)
mutations_data = mutations_data[which(mutations_data$PATIENT_ID%in%clinical_data$PATIENT_ID),]
rownames(mutations_data) = 1:nrow(mutations_data)

# perform the survival association analysis
set.seed(12345)
model_estimates = list()
input_data = list()
for(i in sort(unique(clinical_data$CANCER_TYPE))) {
    # process the data for the current cancer type
    survival_data = clinical_data[which(clinical_data$CANCER_TYPE==i),]
    rownames(survival_data) = survival_data$PATIENT_ID
    curr_mutations_data = mutations_data[which(mutations_data$PATIENT_ID%in%survival_data$PATIENT_ID),]
    mutations = matrix(0,nrow=nrow(survival_data),ncol=length(unique(curr_mutations_data$GENE_NAME)))
    rownames(mutations) = rownames(survival_data)
    colnames(mutations) = sort(unique(curr_mutations_data$GENE_NAME))
    for(j in 1:nrow(curr_mutations_data)) {
        mutations[curr_mutations_data[j,"PATIENT_ID"],curr_mutations_data[j,"GENE_NAME"]] = curr_mutations_data[j,"VAF"]
    }
    genes_counts = colSums(apply(X=mutations,MARGIN=2,FUN=function(x){x>0.00}))
    valid_genes1 = names(which((genes_counts/nrow(mutations))>0.005))
    valid_genes2 = names(which(genes_counts>3))
    features_data = mutations[,sort(unique(intersect(valid_genes1,valid_genes2))),drop=FALSE]
    features_data_mutated = features_data
    features_data_mutated[which(features_data_mutated>0)] = 1
    colnames(features_data) = paste0(colnames(features_data),"_VAF")
    colnames(features_data_mutated) = paste0(colnames(features_data_mutated),"_MUTATED")
    analysis_data = data.frame(cbind(survival_data[,c("OS_MONTHS","OS_STATUS","AGE","GENDER")],features_data,features_data_mutated))
    colnames(analysis_data)[1:2] = c("Times","Status")
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
        NUM_MUT = rep(NA,nrow(model_estimate))
        for(j in 1:nrow(model_estimate)) {
            if(model_estimate$VARIABLE[j]!="AGE"&&model_estimate$VARIABLE[j]!="GENDER") {
                NUM_MUT[j] = sum(features_data_mutated[,gsub("_VAF","_MUTATED",model_estimate$VARIABLE[j])])
            }
        }
        FREQ_MUT = (NUM_MUT/nrow(features_data_mutated))
        model_estimate$NUM_MUT = NUM_MUT
        model_estimate$FREQ_MUT = FREQ_MUT
        model_estimates[[i]] = model_estimate
    }
}
selected_models = NULL
for(i in 1:length(model_estimates)) {
    selected_models = c(selected_models,length(grep("_VAF",model_estimates[[i]][,"VARIABLE"])))
}
model_estimates = model_estimates[which(selected_models>0)]
names(model_estimates)[2] = "Breast Ductal Carcinoma"
names(model_estimates)[7] = "High-Grade Ovarian Cancer"

# save the results
save(input_data,file="results/input_data_vaf.RData")
save(model_estimates,file="results/model_estimates_vaf.RData")
write.xlsx(x=model_estimates,file="results/model_estimates_vaf.xlsx",rowNames=FALSE,colNames=TRUE)
