# load the required libraries
library("data.table")

# process the clinical data
data_clinical_patient = fread(file="downloads/heme_msk_impact_2022/data_clinical_patient.txt",skip=4)
data_clinical_sample = fread(file="downloads/heme_msk_impact_2022/data_clinical_sample.txt",skip=4)
PATIENT_ID = as.character(data_clinical_sample$PATIENT_ID)
SAMPLE_ID = as.character(data_clinical_sample$SAMPLE_ID)
SAMPLE_TYPE = as.character(data_clinical_sample$SAMPLE_TYPE)
SAMPLE_CLASS = as.character(data_clinical_sample$SAMPLE_CLASS)
PRIMARY_SITE = as.character(data_clinical_sample$PRIMARY_SITE)
METASTATIC_SITE = as.character(data_clinical_sample$METASTATIC_SITE)
CANCER_TYPE = as.character(data_clinical_sample$CANCER_TYPE)
CANCER_TYPE_DETAILED = as.character(data_clinical_sample$CANCER_TYPE_DETAILED)
cont = 0
AGE = rep(NA,length(PATIENT_ID))
GENDER = rep(NA,length(PATIENT_ID))
RACE = rep(NA,length(PATIENT_ID))
ETHNICITY = rep(NA,length(PATIENT_ID))
OS_MONTHS = rep(NA,length(PATIENT_ID))
OS_STATUS = rep(NA,length(PATIENT_ID))
for(i in PATIENT_ID) {
    pos = which(data_clinical_patient$PATIENT_ID==i)
    cont = cont + 1
    if(length(pos)==1) {
        GENDER[cont] = as.character(data_clinical_patient$SEX[pos])
        RACE[cont] = as.character(data_clinical_patient$RACE[pos])
        ETHNICITY[cont] = as.character(data_clinical_patient$ETHNICITY[pos])
        OS_MONTHS[cont] = as.numeric(data_clinical_patient$OS_MONTHS[pos])
        OS_STATUS[cont] = as.character(data_clinical_patient$OS_STATUS[pos])
    }
}
clinical_data = data.frame(PATIENT_ID=PATIENT_ID,SAMPLE_ID=SAMPLE_ID,AGE=AGE,GENDER=GENDER,RACE=RACE,ETHNICITY=ETHNICITY,SAMPLE_TYPE=SAMPLE_TYPE,SAMPLE_CLASS=SAMPLE_CLASS,PRIMARY_SITE=PRIMARY_SITE,METASTATIC_SITE=METASTATIC_SITE,CANCER_TYPE=CANCER_TYPE,CANCER_TYPE_DETAILED=CANCER_TYPE_DETAILED,OS_MONTHS=OS_MONTHS,OS_STATUS=OS_STATUS,check.names=FALSE,stringsAsFactors=FALSE)
clinical_data$GENDER[which(clinical_data$GENDER=="")] = NA
clinical_data$RACE[which(clinical_data$RACE%in%c("","NO VALUE ENTERED","PT REFUSED TO ANSWER","UNKNOWN"))] = NA
clinical_data$ETHNICITY[which(clinical_data$ETHNICITY=="")] = NA
clinical_data$SAMPLE_TYPE[which(clinical_data$SAMPLE_TYPE%in%c("","Unknown"))] = NA
clinical_data$SAMPLE_CLASS[which(clinical_data$SAMPLE_CLASS=="")] = NA
clinical_data$PRIMARY_SITE[which(clinical_data$PRIMARY_SITE%in%c("","Unknown"))] = NA
clinical_data$METASTATIC_SITE[which(clinical_data$METASTATIC_SITE%in%c("","Unknown"))] = NA
clinical_data$OS_STATUS[which(clinical_data$OS_STATUS=="")] = NA
clinical_data = clinical_data[order(clinical_data$PATIENT_ID,clinical_data$SAMPLE_ID),]
clinical_data = unique(clinical_data)
curr_duplicated = which(duplicated(clinical_data$PATIENT_ID,fromLast=TRUE))
if(length(curr_duplicated)>0) {
    clinical_data = clinical_data[-curr_duplicated,]
}
rownames(clinical_data) = 1:nrow(clinical_data)

# process the mutations data
data_mutations_extended = fread(file="downloads/heme_msk_impact_2022/data_mutations.txt")
data_mutations_extended = data_mutations_extended[,c("Tumor_Sample_Barcode","Hugo_Symbol","Chromosome","Start_Position","End_Position","Consequence","Variant_Classification","Variant_Type","Reference_Allele","Tumor_Seq_Allele2","t_alt_count","t_ref_count")]
data_mutations_extended = data_mutations_extended[which(data_mutations_extended$Hugo_Symbol!="Unknown"),]
data_mutations_extended$Consequence[which(data_mutations_extended$Consequence=="")] = NA
data_mutations_extended$Variant_Classification[which(data_mutations_extended$Variant_Classification=="")] = NA
data_mutations_extended$Variant_Type[which(data_mutations_extended$Variant_Type=="")] = NA
data_mutations_extended = data_mutations_extended[-which(data_mutations_extended$Variant_Classification=="Intron"),] # remove Intronic mutations
data_mutations_extended = data_mutations_extended[-which(data_mutations_extended$Variant_Classification=="Silent"),] # remove Silent mutations
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
mutations_data = data.frame(SAMPLE_ID=SAMPLE_ID,GENE_NAME=GENE_NAME,CHROMOSOME=CHROMOSOME,START_POSITION=START_POSITION,END_POSITION=END_POSITION,CONSEQUENCE=CONSEQUENCE,VARIANT_CLASSIFICATION=VARIANT_CLASSIFICATION,VARIANT_TYPE=VARIANT_TYPE,REF_ALLELE=REF_ALLELE,ALT_ALLELE=ALT_ALLELE,ALT_COUNT=ALT_COUNT,REF_COUNT=REF_COUNT,check.names=FALSE,stringsAsFactors=FALSE)
mutations_data = mutations_data[order(mutations_data$SAMPLE_ID,mutations_data$GENE_NAME,mutations_data$REF_COUNT,mutations_data$ALT_COUNT),]
mutations_data = unique(mutations_data)
rownames(mutations_data) = 1:nrow(mutations_data)
mutations = matrix(0, nrow = length(unique(mutations_data$SAMPLE_ID)), ncol = length(unique(mutations_data$GENE_NAME)))
rownames(mutations) = sort(unique(mutations_data$SAMPLE_ID))
colnames(mutations) = sort(unique(mutations_data$GENE_NAME))
for(i in rownames(mutations)) {
    mutations[i,mutations_data$GENE_NAME[which(mutations_data$SAMPLE_ID==i)]] = 1
}
mutations = mutations[sort(unique(rownames(mutations))),]
mutations = mutations[,sort(unique(colnames(mutations)))]
mutations = as.matrix(mutations)

# process the copy number data
data_cna_hg19 = fread(file="downloads/heme_msk_impact_2022/data_cna_hg19.seg")
SAMPLE_ID = as.character(data_cna_hg19$ID)
CHROMOSOME = as.character(data_cna_hg19$chrom)
START_POSITION = as.numeric(data_cna_hg19$loc.start)
END_POSITION = as.numeric(data_cna_hg19$loc.end)
SEG_MEAN = as.numeric(data_cna_hg19$seg.mean)
cna_data = data.frame(SAMPLE_ID=SAMPLE_ID,CHROMOSOME=CHROMOSOME,START_POSITION=START_POSITION,END_POSITION=END_POSITION,SEG_MEAN=SEG_MEAN,check.names=FALSE,stringsAsFactors=FALSE)
cna_data = cna_data[order(cna_data$SAMPLE_ID,cna_data$CHROMOSOME,cna_data$START_POSITION,cna_data$END_POSITION,cna_data$SEG_MEAN),]
cna_data = unique(cna_data)
rownames(cna_data) = 1:nrow(cna_data)
data_cna = fread(file="downloads/heme_msk_impact_2022/data_cna.txt")
genes = data_cna$Hugo_Symbol
data_cna$Hugo_Symbol = NULL
gistic = as.matrix(data_cna)
gistic = t(gistic)
colnames(gistic) = genes
gistic = gistic[sort(unique(rownames(gistic))),]
gistic = gistic[,sort(unique(colnames(gistic)))]

# save the results
save(clinical_data,file="processed_data/clinical_data.RData")
save(mutations_data,file="processed_data/mutations_data.RData")
save(mutations,file="processed_data/mutations.RData")
save(cna_data,file="processed_data/cna_data.RData")
save(gistic,file="processed_data/gistic.RData")
