# load the processed data
clinical_data = NULL
mutations_data = NULL
mutations = NULL
load(file="processed_data/2020/clinical_data.RData")
load(file="processed_data/2020/mutations_data.RData")
load(file="processed_data/2020/mutations.RData")
clinical_data_2020 = clinical_data
mutations_data_2020 = mutations_data
mutations_2020 = mutations
clinical_data = NULL
mutations_data = NULL
mutations = NULL
load(file="processed_data/2023/clinical_data.RData")
load(file="processed_data/2023/mutations_data.RData")
load(file="processed_data/2023/mutations.RData")
clinical_data_2023 = clinical_data
mutations_data_2023 = mutations_data
mutations_2023 = mutations
clinical_data = NULL
mutations_data = NULL
mutations = NULL

# select the unique samples
valid_samples_2020 = sort(unique(clinical_data_2020$PATIENT_ID[which(!clinical_data_2020$PATIENT_ID%in%clinical_data_2023$PATIENT_ID)]))
clinical_data_2020 = clinical_data_2020[which(clinical_data_2020$PATIENT_ID%in%valid_samples_2020),]
rownames(clinical_data_2020) = 1:nrow(clinical_data_2020)
mutations_data_2020 = mutations_data_2020[which(mutations_data_2020$SAMPLE_ID%in%clinical_data_2020$SAMPLE_ID),]
rownames(mutations_data_2020) = 1:nrow(mutations_data_2020)
mutations_2020 = mutations_2020[sort(unique(mutations_data_2020$SAMPLE_ID)),]
mutations_2020 = mutations_2020[,sort(unique(names(which(colSums(mutations_2020)>0))))]

# renaming from SAMPLE_IDs to PATIENT_IDs
for(i in 1:nrow(mutations_data_2020)) {
    mutations_data_2020$SAMPLE_ID[i] = clinical_data_2020$PATIENT_ID[which(clinical_data_2020$SAMPLE_ID==mutations_data_2020$SAMPLE_ID[i])]
}
colnames(mutations_data_2020)[1] = "PATIENT_ID"
for(i in 1:nrow(mutations_2020)) {
    rownames(mutations_2020)[i] = clinical_data_2020$PATIENT_ID[which(clinical_data_2020$SAMPLE_ID==rownames(mutations_2020)[i])]
}
for(i in 1:nrow(mutations_data_2023)) {
    mutations_data_2023$SAMPLE_ID[i] = clinical_data_2023$PATIENT_ID[which(clinical_data_2023$SAMPLE_ID==mutations_data_2023$SAMPLE_ID[i])]
}
colnames(mutations_data_2023)[1] = "PATIENT_ID"
for(i in 1:nrow(mutations_2023)) {
    rownames(mutations_2023)[i] = clinical_data_2023$PATIENT_ID[which(clinical_data_2023$SAMPLE_ID==rownames(mutations_2023)[i])]
}

# save the results
save(clinical_data_2020,file="final_data/clinical_data_2020.RData")
save(mutations_data_2020,file="final_data/mutations_data_2020.RData")
save(mutations_2020,file="final_data/mutations_2020.RData")
save(clinical_data_2023,file="final_data/clinical_data_2023.RData")
save(mutations_data_2023,file="final_data/mutations_data_2023.RData")
save(mutations_2023,file="final_data/mutations_2023.RData")
