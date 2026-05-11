# load and process the data
load("results/mutations_data.RData")
muts = paste0(mutations_data$GENE_NAME,"_",mutations_data$HGVSp)
muts = table(muts)
muts = muts[which(muts>=10)]
muts = sort(muts,decreasing=TRUE)
print(muts[1:20])
