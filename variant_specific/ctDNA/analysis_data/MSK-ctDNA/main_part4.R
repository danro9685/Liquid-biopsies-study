# load the required libraries
library("data.table")

# process the clinical data
data_clinical_patient = fread(file="downloads/2022/nsclc_ctdx_msk_2022/data_clinical_patient.txt",skip=4)
data_clinical_sample = fread(file="downloads/2022/nsclc_ctdx_msk_2022/data_clinical_sample.txt",skip=4)
PATIENT_ID = as.character(data_clinical_sample$PATIENT_ID)
SAMPLE_ID = as.character(data_clinical_sample$SAMPLE_ID)
SAMPLE_TYPE = as.character(data_clinical_sample$SAMPLE_TYPE)
CANCER_TYPE = as.character(data_clinical_sample$CANCER_TYPE)
CANCER_TYPE_DETAILED = as.character(data_clinical_sample$CANCER_TYPE_DETAILED)
SAMPLE_CLASS = as.character(data_clinical_sample$SAMPLE_CLASS)
cont = 0
AGE = rep(NA,length(PATIENT_ID))
GENDER = rep(NA,length(PATIENT_ID))
RACE = rep(NA,length(PATIENT_ID))
OS_MONTHS = rep(NA,length(PATIENT_ID))
OS_STATUS = rep(NA,length(PATIENT_ID))
for(i in PATIENT_ID) {
    pos = which(data_clinical_patient$PATIENT_ID==i)
    cont = cont + 1
    if(length(pos)==1) {
        AGE[cont] = as.numeric(data_clinical_patient$AGE_CURRENT[pos])
        GENDER[cont] = as.character(data_clinical_patient$SEX[pos])
        RACE[cont] = as.character(data_clinical_patient$RACE[pos])
        OS_MONTHS[cont] = as.numeric(data_clinical_patient$OS_MONTHS[pos])
        OS_STATUS[cont] = as.character(data_clinical_patient$OS_STATUS[pos])
    }
}
clinical_data = data.frame(PATIENT_ID=PATIENT_ID,SAMPLE_ID=SAMPLE_ID,AGE=AGE,GENDER=GENDER,RACE=RACE,SAMPLE_TYPE=SAMPLE_TYPE,CANCER_TYPE=CANCER_TYPE,CANCER_TYPE_DETAILED=CANCER_TYPE_DETAILED,OS_MONTHS=OS_MONTHS,OS_STATUS=OS_STATUS,SAMPLE_CLASS=SAMPLE_CLASS,check.names=FALSE,stringsAsFactors=FALSE)
clinical_data$GENDER[which(clinical_data$GENDER=="")] = NA
clinical_data$RACE[which(clinical_data$RACE%in%c("","NO VALUE ENTERED","PT REFUSED TO ANSWER","UNKNOWN"))] = NA
clinical_data$OS_STATUS[which(clinical_data$OS_STATUS=="")] = NA
clinical_data = clinical_data[which(clinical_data$SAMPLE_CLASS=="cfDNA"),]
clinical_data$SAMPLE_CLASS = NULL
curr_duplicated = which(duplicated(clinical_data$PATIENT_ID,fromLast=TRUE))
if(length(curr_duplicated)>0) {
    clinical_data = clinical_data[-curr_duplicated,]
}
clinical_data = clinical_data[order(clinical_data$PATIENT_ID,clinical_data$SAMPLE_ID),]
clinical_data = unique(clinical_data)
rownames(clinical_data) = 1:nrow(clinical_data)

# process the mutations data
data_mutations_extended = fread(file="downloads/2022/nsclc_ctdx_msk_2022/data_mutations.txt")
data_mutations_extended = data_mutations_extended[,c("Tumor_Sample_Barcode","Hugo_Symbol","Chromosome","Start_Position","End_Position","Consequence","Variant_Classification","Variant_Type","Reference_Allele","Tumor_Seq_Allele2","t_alt_count","t_ref_count")]
data_mutations_extended = data_mutations_extended[which(data_mutations_extended$Hugo_Symbol!="Unknown"),]
data_mutations_extended$Consequence[which(data_mutations_extended$Consequence=="")] = NA
data_mutations_extended$Variant_Classification[which(data_mutations_extended$Variant_Classification=="")] = NA
data_mutations_extended$Variant_Type[which(data_mutations_extended$Variant_Type=="")] = NA
data_mutations_extended = data_mutations_extended[-which(data_mutations_extended$Variant_Classification=="Intron"),] # remove Intronic mutations
data_mutations_extended = data_mutations_extended[-which(data_mutations_extended$Variant_Classification=="Silent"),] # remove Silent mutations
data_mutations_extended = data_mutations_extended[which(data_mutations_extended$Tumor_Sample_Barcode%in%clinical_data$SAMPLE_ID),]
SAMPLE_ID = as.character(data_mutations_extended$Tumor_Sample_Barcode)
GENE_NAME = as.character(data_mutations_extended$Hugo_Symbol)
CHROMOSOME = as.character(data_mutations_extended$Chromosome)
START_POSITION = as.numeric(data_mutations_extended$Start_Position)
END_POSITION = as.numeric(data_mutations_extended$End_Position)
CONSEQUENCE = as.character(data_mutations_extended$Consequence)
VARIANT_CLASSIFICATION = as.character(data_mutations_extended$Variant_Classification)
VARIANT_TYPE = as.character(data_mutations_extended$Variant_Type)
REF_ALLELE = as.character(data_mutations_extended$Reference_Allele)
ALT_ALLELE = as.character(data_mutations_extended$Tumor_Seq_Allele2)
ALT_COUNT = as.numeric(data_mutations_extended$t_alt_count)
REF_COUNT = as.numeric(data_mutations_extended$t_ref_count)
GENE_NAME[which(GENE_NAME=="CDKN2Ap14ARF")] = "CDKN2A"
GENE_NAME[which(GENE_NAME=="CDKN2Ap16INK4A")] = "CDKN2A"
GENE_NAME[which(GENE_NAME=="MLL")] = "KMT2A"
GENE_NAME[which(GENE_NAME=="MLL2")] = "KMT2D"
GENE_NAME[which(GENE_NAME=="MLL3")] = "KMT2C"
GENE_NAME[which(GENE_NAME=="MLL4")] = "KMT2B"
GENE_NAME[which(GENE_NAME=="SETD8")] = "KMT5A"
mutations_data = data.frame(SAMPLE_ID=SAMPLE_ID,GENE_NAME=GENE_NAME,CHROMOSOME=CHROMOSOME,START_POSITION=START_POSITION,END_POSITION=END_POSITION,CONSEQUENCE=CONSEQUENCE,VARIANT_CLASSIFICATION=VARIANT_CLASSIFICATION,VARIANT_TYPE=VARIANT_TYPE,REF_ALLELE=REF_ALLELE,ALT_ALLELE=ALT_ALLELE,ALT_COUNT=ALT_COUNT,REF_COUNT=REF_COUNT,check.names=FALSE,stringsAsFactors=FALSE)
mutations_data = mutations_data[order(mutations_data$SAMPLE_ID,mutations_data$GENE_NAME,mutations_data$REF_COUNT,mutations_data$ALT_COUNT),]
mutations_data = unique(mutations_data)
rownames(mutations_data) = 1:nrow(mutations_data)
samples = sort(unique(mutations_data$SAMPLE_ID))
genes = sort(unique(mutations_data$GENE_NAME))
mutations = matrix(0, nrow = length(samples), ncol = length(genes))
rownames(mutations) = samples
colnames(mutations) = genes
for(i in 1:nrow(mutations_data)) {
    mutations[mutations_data$SAMPLE_ID[i],mutations_data$GENE_NAME[i]] = 1
}

# process the svs data
data_sv = fread(file="downloads/2022/nsclc_ctdx_msk_2022/data_sv.txt")
SAMPLE_ID = as.character(data_sv$Sample_Id)
GENE_NAME1 = as.character(data_sv$Site1_Hugo_Symbol)
GENE_NAME2 = as.character(data_sv$Site2_Hugo_Symbol)
EVENT_TYPE = as.character(data_sv$Event_Info)
sv_data = data.frame(SAMPLE_ID=SAMPLE_ID,GENE_NAME1=GENE_NAME1,GENE_NAME2=GENE_NAME2,EVENT_TYPE=EVENT_TYPE,check.names=FALSE,stringsAsFactors=FALSE)
sv_data$GENE_NAME2[which(sv_data$GENE_NAME2=="")] = NA
sv_data = unique(sv_data)
sv_data = sv_data[order(sv_data$SAMPLE_ID,sv_data$GENE_NAME1,sv_data$GENE_NAME2),]
sv_data = sv_data[which(sv_data$SAMPLE_ID%in%clinical_data$SAMPLE_ID),]
rownames(sv_data) = 1:nrow(sv_data)

# save the results
save(clinical_data,file="processed_data/2022/cfDNA/clinical_data.RData")
save(mutations_data,file="processed_data/2022/cfDNA/mutations_data.RData")
save(mutations,file="processed_data/2022/cfDNA/mutations.RData")
save(sv_data,file="processed_data/2022/cfDNA/sv_data.RData")
