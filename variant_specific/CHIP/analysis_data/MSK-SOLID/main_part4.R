# load the processed data
clinical_data = NULL
mutations_data = NULL
mutations = NULL
gistic = NULL
load(file="processed_data/2017/clinical_data.RData")
load(file="processed_data/2017/mutations_data.RData")
load(file="processed_data/2017/mutations.RData")
load(file="processed_data/2017/cna_data.RData")
load(file="processed_data/2017/gistic.RData")
clinical_data_2017 = clinical_data
mutations_data_2017 = mutations_data
mutations_2017 = mutations
cna_data_2017 = cna_data
gistic_2017 = gistic
clinical_data = NULL
mutations_data = NULL
mutations = NULL
cna_data = NULL
gistic = NULL
load(file="processed_data/2021/clinical_data.RData")
load(file="processed_data/2021/mutations_data.RData")
load(file="processed_data/2021/mutations.RData")
load(file="processed_data/2021/cna_data.RData")
load(file="processed_data/2021/gistic.RData")
clinical_data_2021 = clinical_data
mutations_data_2021 = mutations_data
mutations_2021 = mutations
cna_data_2021 = cna_data
gistic_2021 = gistic
clinical_data = NULL
mutations_data = NULL
mutations = NULL
cna_data = NULL
gistic = NULL
load(file="processed_data/2024/clinical_data.RData")
load(file="processed_data/2024/mutations_data.RData")
load(file="processed_data/2024/mutations.RData")
load(file="processed_data/2024/cna_data.RData")
load(file="processed_data/2024/gistic.RData")
clinical_data_2024 = clinical_data
mutations_data_2024 = mutations_data
mutations_2024 = mutations
cna_data_2024 = cna_data
gistic_2024 = gistic
clinical_data = NULL
mutations_data = NULL
mutations = NULL
cna_data = NULL
gistic = NULL

# select the samples with all the omics
valid_samples = sort(unique(intersect(clinical_data_2017$SAMPLE_ID,mutations_data_2017$SAMPLE_ID)))
valid_samples = sort(unique(intersect(valid_samples,rownames(mutations_2017))))
valid_samples = sort(unique(intersect(valid_samples,cna_data_2017$SAMPLE_ID)))
valid_samples = sort(unique(intersect(valid_samples,rownames(gistic_2017))))
clinical_data_2017 = clinical_data_2017[which(clinical_data_2017$SAMPLE_ID%in%valid_samples),]
rownames(clinical_data_2017) = 1:nrow(clinical_data_2017)
mutations_data_2017 = mutations_data_2017[which(mutations_data_2017$SAMPLE_ID%in%clinical_data_2017$SAMPLE_ID),]
rownames(mutations_data_2017) = 1:nrow(mutations_data_2017)
mutations_2017 = mutations_2017[clinical_data_2017$SAMPLE_ID,]
cna_data_2017 = cna_data_2017[which(cna_data_2017$SAMPLE_ID%in%valid_samples),]
rownames(cna_data_2017) = 1:nrow(cna_data_2017)
gistic_2017 = gistic_2017[clinical_data_2017$SAMPLE_ID,]
valid_samples = sort(unique(intersect(clinical_data_2021$SAMPLE_ID,mutations_data_2021$SAMPLE_ID)))
valid_samples = sort(unique(intersect(valid_samples,rownames(mutations_2021))))
valid_samples = sort(unique(intersect(valid_samples,cna_data_2021$SAMPLE_ID)))
valid_samples = sort(unique(intersect(valid_samples,rownames(gistic_2021))))
clinical_data_2021 = clinical_data_2021[which(clinical_data_2021$SAMPLE_ID%in%valid_samples),]
rownames(clinical_data_2021) = 1:nrow(clinical_data_2021)
mutations_data_2021 = mutations_data_2021[which(mutations_data_2021$SAMPLE_ID%in%clinical_data_2021$SAMPLE_ID),]
rownames(mutations_data_2021) = 1:nrow(mutations_data_2021)
mutations_2021 = mutations_2021[clinical_data_2021$SAMPLE_ID,]
cna_data_2021 = cna_data_2021[which(cna_data_2021$SAMPLE_ID%in%valid_samples),]
rownames(cna_data_2021) = 1:nrow(cna_data_2021)
gistic_2021 = gistic_2021[clinical_data_2021$SAMPLE_ID,]
valid_samples = sort(unique(intersect(clinical_data_2024$SAMPLE_ID,mutations_data_2024$SAMPLE_ID)))
valid_samples = sort(unique(intersect(valid_samples,rownames(mutations_2024))))
valid_samples = sort(unique(intersect(valid_samples,cna_data_2024$SAMPLE_ID)))
valid_samples = sort(unique(intersect(valid_samples,rownames(gistic_2024))))
clinical_data_2024 = clinical_data_2024[which(clinical_data_2024$SAMPLE_ID%in%valid_samples),]
rownames(clinical_data_2024) = 1:nrow(clinical_data_2024)
mutations_data_2024 = mutations_data_2024[which(mutations_data_2024$SAMPLE_ID%in%clinical_data_2024$SAMPLE_ID),]
rownames(mutations_data_2024) = 1:nrow(mutations_data_2024)
mutations_2024 = mutations_2024[clinical_data_2024$SAMPLE_ID,]
cna_data_2024 = cna_data_2024[which(cna_data_2024$SAMPLE_ID%in%valid_samples),]
rownames(cna_data_2024) = 1:nrow(cna_data_2024)
gistic_2024 = gistic_2024[clinical_data_2024$SAMPLE_ID,]

# select the unique samples
valid_samples_2017 = sort(unique(clinical_data_2017$PATIENT_ID[which(!clinical_data_2017$PATIENT_ID%in%clinical_data_2021$PATIENT_ID)]))
valid_samples_2017 = sort(unique(valid_samples_2017[which(!valid_samples_2017%in%clinical_data_2024$PATIENT_ID)]))
clinical_data_2017 = clinical_data_2017[which(clinical_data_2017$PATIENT_ID%in%valid_samples_2017),]
rownames(clinical_data_2017) = 1:nrow(clinical_data_2017)
mutations_data_2017 = mutations_data_2017[which(mutations_data_2017$SAMPLE_ID%in%clinical_data_2017$SAMPLE_ID),]
rownames(mutations_data_2017) = 1:nrow(mutations_data_2017)
mutations_2017 = mutations_2017[sort(unique(mutations_data_2017$SAMPLE_ID)),]
mutations_2017 = mutations_2017[,sort(unique(names(which(colSums(mutations_2017)>0))))]
cna_data_2017 = cna_data_2017[which(cna_data_2017$SAMPLE_ID%in%clinical_data_2017$SAMPLE_ID),]
rownames(cna_data_2017) = 1:nrow(cna_data_2017)
gistic_2017 = gistic_2017[sort(unique(clinical_data_2017$SAMPLE_ID)),]
valid_samples_2021 = sort(unique(clinical_data_2021$PATIENT_ID[which(!clinical_data_2021$PATIENT_ID%in%clinical_data_2024$PATIENT_ID)]))
clinical_data_2021 = clinical_data_2021[which(clinical_data_2021$PATIENT_ID%in%valid_samples_2021),]
rownames(clinical_data_2021) = 1:nrow(clinical_data_2021)
mutations_data_2021 = mutations_data_2021[which(mutations_data_2021$SAMPLE_ID%in%clinical_data_2021$SAMPLE_ID),]
rownames(mutations_data_2021) = 1:nrow(mutations_data_2021)
mutations_2021 = mutations_2021[sort(unique(mutations_data_2021$SAMPLE_ID)),]
mutations_2021 = mutations_2021[,sort(unique(names(which(colSums(mutations_2021)>0))))]
cna_data_2021 = cna_data_2021[which(cna_data_2021$SAMPLE_ID%in%clinical_data_2021$SAMPLE_ID),]
rownames(cna_data_2021) = 1:nrow(cna_data_2021)
gistic_2021 = gistic_2021[sort(unique(clinical_data_2021$SAMPLE_ID)),]

# remove missing data
clinical_data_2024$SUBTYPE = NA
clinical_data_2024 = clinical_data_2024[,colnames(clinical_data_2021)]
is.missing = (colSums(apply(X=gistic_2024,MARGIN=2,FUN=is.na))/nrow(gistic_2024))
gistic_2024 = gistic_2024[,sort(unique(names(which(is.missing<0.01))))]
gistic_2024[which(is.na(gistic_2024))] = 0

# renaming from SAMPLE_IDs to PATIENT_IDs
for(i in 1:nrow(mutations_data_2017)) {
    mutations_data_2017$SAMPLE_ID[i] = clinical_data_2017$PATIENT_ID[which(clinical_data_2017$SAMPLE_ID==mutations_data_2017$SAMPLE_ID[i])]
}
colnames(mutations_data_2017)[1] = "PATIENT_ID"
for(i in 1:nrow(mutations_2017)) {
    rownames(mutations_2017)[i] = clinical_data_2017$PATIENT_ID[which(clinical_data_2017$SAMPLE_ID==rownames(mutations_2017)[i])]
}
for(i in 1:nrow(cna_data_2017)) {
    cna_data_2017$SAMPLE_ID[i] = clinical_data_2017$PATIENT_ID[which(clinical_data_2017$SAMPLE_ID==cna_data_2017$SAMPLE_ID[i])]
}
colnames(cna_data_2017)[1] = "PATIENT_ID"
for(i in 1:nrow(gistic_2017)) {
    rownames(gistic_2017)[i] = clinical_data_2017$PATIENT_ID[which(clinical_data_2017$SAMPLE_ID==rownames(gistic_2017)[i])]
}
for(i in 1:nrow(mutations_data_2021)) {
    mutations_data_2021$SAMPLE_ID[i] = clinical_data_2021$PATIENT_ID[which(clinical_data_2021$SAMPLE_ID==mutations_data_2021$SAMPLE_ID[i])]
}
colnames(mutations_data_2021)[1] = "PATIENT_ID"
for(i in 1:nrow(mutations_2021)) {
    rownames(mutations_2021)[i] = clinical_data_2021$PATIENT_ID[which(clinical_data_2021$SAMPLE_ID==rownames(mutations_2021)[i])]
}
for(i in 1:nrow(cna_data_2021)) {
    cna_data_2021$SAMPLE_ID[i] = clinical_data_2021$PATIENT_ID[which(clinical_data_2021$SAMPLE_ID==cna_data_2021$SAMPLE_ID[i])]
}
colnames(cna_data_2021)[1] = "PATIENT_ID"
for(i in 1:nrow(gistic_2021)) {
    rownames(gistic_2021)[i] = clinical_data_2021$PATIENT_ID[which(clinical_data_2021$SAMPLE_ID==rownames(gistic_2021)[i])]
}
for(i in 1:nrow(mutations_data_2024)) {
    mutations_data_2024$SAMPLE_ID[i] = clinical_data_2024$PATIENT_ID[which(clinical_data_2024$SAMPLE_ID==mutations_data_2024$SAMPLE_ID[i])]
}
colnames(mutations_data_2024)[1] = "PATIENT_ID"
for(i in 1:nrow(mutations_2024)) {
    rownames(mutations_2024)[i] = clinical_data_2024$PATIENT_ID[which(clinical_data_2024$SAMPLE_ID==rownames(mutations_2024)[i])]
}
for(i in 1:nrow(cna_data_2024)) {
    cna_data_2024$SAMPLE_ID[i] = clinical_data_2024$PATIENT_ID[which(clinical_data_2024$SAMPLE_ID==cna_data_2024$SAMPLE_ID[i])]
}
colnames(cna_data_2024)[1] = "PATIENT_ID"
for(i in 1:nrow(gistic_2024)) {
    rownames(gistic_2024)[i] = clinical_data_2024$PATIENT_ID[which(clinical_data_2024$SAMPLE_ID==rownames(gistic_2024)[i])]
}

# save the results
save(clinical_data_2017,file="final_data/clinical_data_2017.RData")
save(mutations_data_2017,file="final_data/mutations_data_2017.RData")
save(mutations_2017,file="final_data/mutations_2017.RData")
save(cna_data_2017,file="final_data/cna_data_2017.RData")
save(gistic_2017,file="final_data/gistic_2017.RData")
save(clinical_data_2021,file="final_data/clinical_data_2021.RData")
save(mutations_data_2021,file="final_data/mutations_data_2021.RData")
save(mutations_2021,file="final_data/mutations_2021.RData")
save(cna_data_2021,file="final_data/cna_data_2021.RData")
save(gistic_2021,file="final_data/gistic_2021.RData")
save(clinical_data_2024,file="final_data/clinical_data_2024.RData")
save(mutations_data_2024,file="final_data/mutations_data_2024.RData")
save(mutations_2024,file="final_data/mutations_2024.RData")
save(cna_data_2024,file="final_data/cna_data_2024.RData")
save(gistic_2024,file="final_data/gistic_2024.RData")
