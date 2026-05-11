# load the required libraries
library("data.table")
library("openxlsx")

# load and process the input data
load("analysis_data/MSK-SOLID/final_data/clinical_data_2021.RData")
load("analysis_data/MSK-SOLID/final_data/clinical_data_2024.RData")
clinical_data = rbind(clinical_data_2021,clinical_data_2024)
clinical_data = clinical_data[order(clinical_data$PATIENT_ID),]
clinical_data_2021 = NULL
clinical_data_2024 = NULL
load("analysis_data/MSK-ctDNA/final_data/mutations_data_2024.RData")
mutations_data = mutations_data_2024
mutations_data_2024 = NULL
valid_samples = sort(unique(intersect(clinical_data$PATIENT_ID,mutations_data$PATIENT_ID)))
clinical_data = clinical_data[which(clinical_data$PATIENT_ID%in%valid_samples),]
rownames(clinical_data) = 1:nrow(clinical_data)
mutations_data = mutations_data[which(mutations_data$PATIENT_ID%in%clinical_data$PATIENT_ID),]
VAF = (mutations_data$ALT_COUNT/(mutations_data$ALT_COUNT+mutations_data$REF_COUNT))
mutations_data$VAF = VAF
mutations_data = mutations_data[which(mutations_data$VAF>0),]
rownames(mutations_data) = 1:nrow(mutations_data)
CANCER_TYPE = rep(NA,nrow(mutations_data))
for(i in 1:nrow(clinical_data)) {
    CANCER_TYPE[which(mutations_data$PATIENT_ID==clinical_data$PATIENT_ID[i])] = clinical_data$CANCER_TYPE_DETAILED[i]
}
mutations_data$CANCER_TYPE = CANCER_TYPE
mutations_data$CANCER_TYPE[which(mutations_data$CANCER_TYPE=="Uterine Carcinosarcoma/Uterine Malignant Mixed Mullerian Tumor")] = "Uterine Carcinosarcoma"
mutations_data$CANCER_TYPE[which(mutations_data$CANCER_TYPE=="Uterine Serous Carcinoma/Uterine Papillary Serous Carcinoma")] = "Uterine Serous Carcinoma"

# process the protein level mutations data
data_mutations_extended = fread(file="analysis_data/MSK-ctDNA/downloads/2024/msk_ctdna_vte_2024/data_mutations.txt")
data_mutations_extended = data_mutations_extended[,c("Tumor_Sample_Barcode","Hugo_Symbol","Chromosome","Start_Position","End_Position","Consequence","Variant_Classification","Variant_Type","Reference_Allele","Tumor_Seq_Allele2","t_alt_count","t_ref_count","HGVSp_Short")]
data_mutations_extended$Tumor_Sample_Barcode = substr(data_mutations_extended$Tumor_Sample_Barcode,1,9)
data_mutations_extended = data_mutations_extended[which(data_mutations_extended$Hugo_Symbol!="Unknown"),]
data_mutations_extended$Consequence[which(data_mutations_extended$Consequence=="")] = NA
data_mutations_extended$Variant_Classification[which(data_mutations_extended$Variant_Classification=="")] = NA
data_mutations_extended$Variant_Type[which(data_mutations_extended$Variant_Type=="")] = NA
data_mutations_extended = data_mutations_extended[-which(data_mutations_extended$Variant_Classification=="Intron"),] # remove Intronic mutations
data_mutations_extended = data_mutations_extended[-which(data_mutations_extended$Variant_Classification=="nonsynonymous_SNV"),] # remove nonsynonymous mutations
data_mutations_extended$Hugo_Symbol[which(data_mutations_extended$Hugo_Symbol=="CDKN2Ap14ARF")] = "CDKN2A"
data_mutations_extended$Hugo_Symbol[which(data_mutations_extended$Hugo_Symbol=="CDKN2Ap16INK4A")] = "CDKN2A"
data_mutations_extended$Hugo_Symbol[which(data_mutations_extended$Hugo_Symbol=="MLL")] = "KMT2A"
data_mutations_extended$Hugo_Symbol[which(data_mutations_extended$Hugo_Symbol=="MLL2")] = "KMT2D"
data_mutations_extended$Hugo_Symbol[which(data_mutations_extended$Hugo_Symbol=="MLL3")] = "KMT2C"
data_mutations_extended$Hugo_Symbol[which(data_mutations_extended$Hugo_Symbol=="MLL4")] = "KMT2B"
data_mutations_extended$Hugo_Symbol[which(data_mutations_extended$Hugo_Symbol=="SETD8")] = "KMT5A"
data_mutations_extended = data_mutations_extended[order(data_mutations_extended$Tumor_Sample_Barcode,data_mutations_extended$Hugo_Symbol,data_mutations_extended$t_ref_count,data_mutations_extended$t_alt_count),]
data_mutations_extended = unique(data_mutations_extended)
rownames(data_mutations_extended) = 1:nrow(data_mutations_extended)
data_mutations_extended = as.data.frame(data_mutations_extended)
ids1 = paste0(mutations_data[,1],mutations_data[,2],mutations_data[,3],mutations_data[,4],mutations_data[,5],mutations_data[,6],mutations_data[,7],mutations_data[,8],mutations_data[,9],mutations_data[,10],mutations_data[,11])
ids2 = paste0(data_mutations_extended[,1],data_mutations_extended[,2],data_mutations_extended[,3],data_mutations_extended[,4],data_mutations_extended[,5],data_mutations_extended[,6],data_mutations_extended[,7],data_mutations_extended[,8],data_mutations_extended[,9],data_mutations_extended[,10],data_mutations_extended[,11])
data_mutations_extended = data_mutations_extended[which(ids2%in%ids1),]
rownames(data_mutations_extended) = 1:nrow(data_mutations_extended)
ids2 = paste0(data_mutations_extended[,1],data_mutations_extended[,2],data_mutations_extended[,3],data_mutations_extended[,4],data_mutations_extended[,5],data_mutations_extended[,6],data_mutations_extended[,7],data_mutations_extended[,8],data_mutations_extended[,9],data_mutations_extended[,10],data_mutations_extended[,11])
print(all(ids1==ids2))
mutations_data$HGVSp = data_mutations_extended$HGVSp_Short

# process the survival data
clinical_data = clinical_data[which(clinical_data$PATIENT_ID%in%mutations_data$PATIENT_ID),]
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
clinical_data$CANCER_TYPE[which(clinical_data$CANCER_TYPE=="Uterine Carcinosarcoma/Uterine Malignant Mixed Mullerian Tumor")] = "Uterine Carcinosarcoma"
clinical_data$CANCER_TYPE[which(clinical_data$CANCER_TYPE=="Uterine Serous Carcinoma/Uterine Papillary Serous Carcinoma")] = "Uterine Serous Carcinoma"

# save the results
save(clinical_data,file="results/clinical_data.RData")
save(mutations_data,file="results/mutations_data.RData")
res = list()
res[["ctDNA mutations"]] = mutations_data
write.xlsx(x=res,file="results/mutations_data_chip.xlsx",rowNames=FALSE,colNames=TRUE)
