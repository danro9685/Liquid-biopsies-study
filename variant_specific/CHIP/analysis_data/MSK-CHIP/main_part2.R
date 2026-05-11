# load the required libraries
library("data.table")

# process the clinical data
data_clinical_sample = fread(file="downloads/2023/msk_ch_2023/data_clinical_sample.txt",skip=4)
PATIENT_ID = as.character(data_clinical_sample$PATIENT_ID)
SAMPLE_ID = as.character(data_clinical_sample$SAMPLE_ID)
CANCER_TYPE = as.character(data_clinical_sample$CANCER_TYPE)
CANCER_TYPE_DETAILED = as.character(data_clinical_sample$CANCER_TYPE_DETAILED)
clinical_data = data.frame(PATIENT_ID=PATIENT_ID,SAMPLE_ID=SAMPLE_ID,CANCER_TYPE=CANCER_TYPE,CANCER_TYPE_DETAILED=CANCER_TYPE_DETAILED,check.names=FALSE,stringsAsFactors=FALSE)
clinical_data = clinical_data[order(clinical_data$PATIENT_ID,clinical_data$SAMPLE_ID),]
clinical_data = unique(clinical_data)
rownames(clinical_data) = 1:nrow(clinical_data)

# process the mutations data
data_mutations_extended = fread(file="downloads/2023/msk_ch_2023/data_mutations.txt")
data_mutations_extended = data_mutations_extended[,c("Tumor_Sample_Barcode","Hugo_Symbol","Chromosome","Start_Position","End_Position","Variant_Classification","Variant_Type","Reference_Allele","Tumor_Seq_Allele2","t_alt_count","t_ref_count")]
data_mutations_extended = data_mutations_extended[which(data_mutations_extended$Hugo_Symbol!="Unknown"),]
data_mutations_extended$Variant_Classification[which(data_mutations_extended$Variant_Classification=="")] = NA
data_mutations_extended$Variant_Type[which(data_mutations_extended$Variant_Type=="")] = NA
data_mutations_extended = data_mutations_extended[-which(data_mutations_extended$Variant_Classification=="Intron"),] # remove Intronic mutations
data_mutations_extended = data_mutations_extended[-which(data_mutations_extended$Variant_Classification=="Silent"),] # remove Silent mutations
SAMPLE_ID = as.character(data_mutations_extended$Tumor_Sample_Barcode)
GENE_NAME = as.character(data_mutations_extended$Hugo_Symbol)
CHROMOSOME = as.character(data_mutations_extended$Chromosome)
START_POSITION = as.numeric(data_mutations_extended$Start_Position)
END_POSITION = as.numeric(data_mutations_extended$End_Position)
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
mutations_data = data.frame(SAMPLE_ID=SAMPLE_ID,GENE_NAME=GENE_NAME,CHROMOSOME=CHROMOSOME,START_POSITION=START_POSITION,END_POSITION=END_POSITION,VARIANT_CLASSIFICATION=VARIANT_CLASSIFICATION,VARIANT_TYPE=VARIANT_TYPE,REF_ALLELE=REF_ALLELE,ALT_ALLELE=ALT_ALLELE,ALT_COUNT=ALT_COUNT,REF_COUNT=REF_COUNT,check.names=FALSE,stringsAsFactors=FALSE)
mutations_data = mutations_data[order(mutations_data$SAMPLE_ID,mutations_data$GENE_NAME,mutations_data$REF_COUNT,mutations_data$ALT_COUNT),]
mutations_data = unique(mutations_data)
rownames(mutations_data) = 1:nrow(mutations_data)

# make the analysis data
valid_samples = intersect(clinical_data$SAMPLE_ID,mutations_data$SAMPLE_ID)
clinical_data = clinical_data[which(clinical_data$SAMPLE_ID%in%valid_samples),]
rownames(clinical_data) = 1:nrow(clinical_data)
mutations_data = mutations_data[which(mutations_data$SAMPLE_ID%in%clinical_data$SAMPLE_ID),]
rownames(mutations_data) = 1:nrow(mutations_data)

# make the binary mutations matrix
samples = sort(unique(mutations_data$SAMPLE_ID))
genes = sort(unique(mutations_data$GENE_NAME))
mutations = matrix(0, nrow = length(samples), ncol = length(genes))
rownames(mutations) = samples
colnames(mutations) = genes
for(i in 1:nrow(mutations_data)) {
    mutations[mutations_data$SAMPLE_ID[i],mutations_data$GENE_NAME[i]] = 1
}

# save the results
save(clinical_data,file="processed_data/2023/clinical_data.RData")
save(mutations_data,file="processed_data/2023/mutations_data.RData")
save(mutations,file="processed_data/2023/mutations.RData")
