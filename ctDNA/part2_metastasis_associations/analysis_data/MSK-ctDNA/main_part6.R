# load the processed data
clinical_data = NULL
mutations_data = NULL
mutations = NULL
load(file="processed_data/2021/cfDNA/clinical_data.RData")
load(file="processed_data/2021/cfDNA/mutations_data.RData")
load(file="processed_data/2021/cfDNA/mutations.RData")
clinical_data_2021 = clinical_data
mutations_data_2021 = mutations_data
mutations_2021 = mutations
clinical_data = NULL
mutations_data = NULL
mutations = NULL
load(file="processed_data/2022/cfDNA/clinical_data.RData")
load(file="processed_data/2022/cfDNA/mutations_data.RData")
load(file="processed_data/2022/cfDNA/mutations.RData")
clinical_data_2022 = clinical_data
mutations_data_2022 = mutations_data
mutations_2022 = mutations
clinical_data = NULL
mutations_data = NULL
mutations = NULL
load(file="processed_data/2024/clinical_data.RData")
load(file="processed_data/2024/mutations_data.RData")
load(file="processed_data/2024/mutations.RData")
clinical_data_2024 = clinical_data
mutations_data_2024 = mutations_data
mutations_2024 = mutations
clinical_data = NULL
mutations_data = NULL
mutations = NULL

# select the samples with complete data
valid_samples = sort(unique(intersect(clinical_data_2021$SAMPLE_ID,mutations_data_2021$SAMPLE_ID)))
valid_samples = sort(unique(intersect(valid_samples,rownames(mutations_2021))))
clinical_data_2021 = clinical_data_2021[which(clinical_data_2021$SAMPLE_ID%in%valid_samples),]
rownames(clinical_data_2021) = 1:nrow(clinical_data_2021)
mutations_data_2021 = mutations_data_2021[which(mutations_data_2021$SAMPLE_ID%in%clinical_data_2021$SAMPLE_ID),]
rownames(mutations_data_2021) = 1:nrow(mutations_data_2021)
mutations_2021 = mutations_2021[clinical_data_2021$SAMPLE_ID,]
valid_samples = sort(unique(intersect(clinical_data_2022$SAMPLE_ID,mutations_data_2022$SAMPLE_ID)))
valid_samples = sort(unique(intersect(valid_samples,rownames(mutations_2022))))
clinical_data_2022 = clinical_data_2022[which(clinical_data_2022$SAMPLE_ID%in%valid_samples),]
rownames(clinical_data_2022) = 1:nrow(clinical_data_2022)
mutations_data_2022 = mutations_data_2022[which(mutations_data_2022$SAMPLE_ID%in%clinical_data_2022$SAMPLE_ID),]
rownames(mutations_data_2022) = 1:nrow(mutations_data_2022)
mutations_2022 = mutations_2022[clinical_data_2022$SAMPLE_ID,]
valid_samples = sort(unique(intersect(clinical_data_2024$SAMPLE_ID,mutations_data_2024$SAMPLE_ID)))
valid_samples = sort(unique(intersect(valid_samples,rownames(mutations_2024))))
clinical_data_2024 = clinical_data_2024[which(clinical_data_2024$SAMPLE_ID%in%valid_samples),]
rownames(clinical_data_2024) = 1:nrow(clinical_data_2024)
clinical_data_2024 = clinical_data_2024[which(clinical_data_2024$SAMPLE_ID%in%clinical_data_2024$SAMPLE_ID),]
rownames(clinical_data_2024) = 1:nrow(clinical_data_2024)
mutations_2024 = mutations_2024[clinical_data_2024$SAMPLE_ID,]

# select the unique samples
valid_samples_2021 = sort(unique(clinical_data_2021$PATIENT_ID[which(!clinical_data_2021$PATIENT_ID%in%clinical_data_2024$PATIENT_ID)]))
valid_samples_2022 = sort(unique(clinical_data_2022$PATIENT_ID[which(!clinical_data_2022$PATIENT_ID%in%clinical_data_2024$PATIENT_ID)]))
valid_samples_2021 = valid_samples_2021[which(!valid_samples_2021%in%valid_samples_2022)]
valid_samples = valid_samples_2021
clinical_data_2021 = clinical_data_2021[which(clinical_data_2021$PATIENT_ID%in%valid_samples),]
rownames(clinical_data_2021) = 1:nrow(clinical_data_2021)
mutations_data_2021 = mutations_data_2021[which(mutations_data_2021$SAMPLE_ID%in%clinical_data_2021$SAMPLE_ID),]
rownames(mutations_data_2021) = 1:nrow(mutations_data_2021)
mutations_2021 = mutations_2021[clinical_data_2021$SAMPLE_ID,]
mutations_2021 = mutations_2021[,sort(unique(names(which(sort(colSums(mutations_2021))>0))))]
valid_samples = valid_samples_2022
clinical_data_2022 = clinical_data_2022[which(clinical_data_2022$PATIENT_ID%in%valid_samples),]
rownames(clinical_data_2022) = 1:nrow(clinical_data_2022)
mutations_data_2022 = mutations_data_2022[which(mutations_data_2022$SAMPLE_ID%in%clinical_data_2022$SAMPLE_ID),]
rownames(mutations_data_2022) = 1:nrow(mutations_data_2022)
mutations_2022 = mutations_2022[clinical_data_2022$SAMPLE_ID,]
mutations_2022 = mutations_2022[,sort(unique(names(which(sort(colSums(mutations_2022))>0))))]
mutations_2024 = mutations_2024[,sort(unique(names(which(sort(colSums(mutations_2024))>0))))]

# renaming from SAMPLE_IDs to PATIENT_IDs
for(i in 1:nrow(mutations_data_2021)) {
    mutations_data_2021$SAMPLE_ID[i] = clinical_data_2021$PATIENT_ID[which(clinical_data_2021$SAMPLE_ID==mutations_data_2021$SAMPLE_ID[i])]
}
colnames(mutations_data_2021)[1] = "PATIENT_ID"
for(i in 1:nrow(mutations_2021)) {
    rownames(mutations_2021)[i] = clinical_data_2021$PATIENT_ID[which(clinical_data_2021$SAMPLE_ID==rownames(mutations_2021)[i])]
}
for(i in 1:nrow(mutations_data_2022)) {
    mutations_data_2022$SAMPLE_ID[i] = clinical_data_2022$PATIENT_ID[which(clinical_data_2022$SAMPLE_ID==mutations_data_2022$SAMPLE_ID[i])]
}
colnames(mutations_data_2022)[1] = "PATIENT_ID"
for(i in 1:nrow(mutations_2022)) {
    rownames(mutations_2022)[i] = clinical_data_2022$PATIENT_ID[which(clinical_data_2022$SAMPLE_ID==rownames(mutations_2022)[i])]
}
for(i in 1:nrow(mutations_data_2024)) {
    mutations_data_2024$SAMPLE_ID[i] = clinical_data_2024$PATIENT_ID[which(clinical_data_2024$SAMPLE_ID==mutations_data_2024$SAMPLE_ID[i])]
}
colnames(mutations_data_2024)[1] = "PATIENT_ID"
for(i in 1:nrow(mutations_2024)) {
    rownames(mutations_2024)[i] = clinical_data_2024$PATIENT_ID[which(clinical_data_2024$SAMPLE_ID==rownames(mutations_2024)[i])]
}

# save the results
save(clinical_data_2021,file="final_data/clinical_data_2021.RData")
save(mutations_data_2021,file="final_data/mutations_data_2021.RData")
save(mutations_2021,file="final_data/mutations_2021.RData")
save(clinical_data_2022,file="final_data/clinical_data_2022.RData")
save(mutations_data_2022,file="final_data/mutations_data_2022.RData")
save(mutations_2022,file="final_data/mutations_2022.RData")
save(clinical_data_2024,file="final_data/clinical_data_2024.RData")
save(mutations_data_2024,file="final_data/mutations_data_2024.RData")
save(mutations_2024,file="final_data/mutations_2024.RData")
