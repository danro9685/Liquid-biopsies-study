# load the required libraries
library("openxlsx")

# load the data for MSK-HEME
load(file="processed_data/clinical_data.RData")
load(file="processed_data/mutations_data.RData")

# load the data for MSK-CHIP
load(file="MSK-CH/clinical_data_2023.RData")
load(file="MSK-CH/mutations_data_2023.RData")

# load the data for MSK-SOLID
load(file="MSK-SOLID/clinical_data_2021.RData")
load(file="MSK-SOLID/mutations_data_2021.RData")
load(file="MSK-SOLID/clinical_data_2024.RData")
load(file="MSK-SOLID/mutations_data_2024.RData")
clinical_data_solid = rbind(clinical_data_2021,clinical_data_2024)
clinical_data_solid = clinical_data_solid[order(clinical_data_solid$PATIENT_ID),]
clinical_data_2021 = NULL
clinical_data_2024 = NULL
mutations_data_solid = rbind(mutations_data_2021,mutations_data_2024)
mutations_data_solid = mutations_data_solid[order(mutations_data_solid$PATIENT_ID),]
mutations_data_2021 = NULL
mutations_data_2024 = NULL

# select the sample with both CHIP and HEME tumors
clinical_data = clinical_data[which(clinical_data$PATIENT_ID%in%intersect(clinical_data$PATIENT_ID,clinical_data_2023$PATIENT_ID)),]
clinical_data = clinical_data[which(clinical_data$SAMPLE_ID%in%intersect(clinical_data$SAMPLE_ID,mutations_data$SAMPLE_ID)),]
rownames(clinical_data) = 1:nrow(clinical_data)
mutations_data = mutations_data[which(mutations_data$SAMPLE_ID%in%clinical_data$SAMPLE_ID),]
rownames(mutations_data) = 1:nrow(mutations_data)

# renaming of the HEME data
clinical_data_heme = clinical_data
mutations_data_heme = mutations_data
for(i in 1:nrow(mutations_data_heme)) {
    mutations_data_heme$SAMPLE_ID[i] = clinical_data_heme$PATIENT_ID[which(clinical_data_heme$SAMPLE_ID==mutations_data_heme$SAMPLE_ID[i])]
}
colnames(mutations_data_heme)[1] = "PATIENT_ID"

# renaming of the CHIP data
clinical_data_ch = clinical_data_2023[which(clinical_data_2023$PATIENT_ID%in%clinical_data_heme$PATIENT_ID),]
rownames(clinical_data_ch) = 1:nrow(clinical_data_ch)
mutations_data_ch = mutations_data_2023[which(mutations_data_2023$PATIENT_ID%in%clinical_data_ch$PATIENT_ID),]
rownames(mutations_data_ch) = 1:nrow(mutations_data_ch)

# renaming of the solid tumors data
clinical_data_solid = clinical_data_solid[which(clinical_data_solid$PATIENT_ID%in%clinical_data_heme$PATIENT_ID),]
rownames(clinical_data_solid) = 1:nrow(clinical_data_solid)
mutations_data_solid = mutations_data_solid[which(mutations_data_solid$PATIENT_ID%in%clinical_data_solid$PATIENT_ID),]
rownames(mutations_data_solid) = 1:nrow(mutations_data_solid)

# save the results
save(clinical_data_ch,file="final_data/clinical_data_ch.RData")
save(mutations_data_ch,file="final_data/mutations_data_ch.RData")
save(clinical_data_heme,file="final_data/clinical_data_heme.RData")
save(mutations_data_heme,file="final_data/mutations_data_heme.RData")
save(clinical_data_solid,file="final_data/clinical_data_solid.RData")
save(mutations_data_solid,file="final_data/mutations_data_solid.RData")
res = list()
res[["Clinical data CHIP"]] = clinical_data_ch
res[["Clinical data HEME"]] = clinical_data_heme
res[["Clinical data solid tumors"]] = clinical_data_solid
res[["Mutations CHIP"]] = mutations_data_ch
res[["Mutations HEME"]] = mutations_data_heme
res[["Mutations solid tumors"]] = mutations_data_solid
write.xlsx(x=res,file="results/Supplementary Table 13 - Clonal dynamics of CHIP in hematologic malignancies.xlsx",rowNames=FALSE,colNames=TRUE)
