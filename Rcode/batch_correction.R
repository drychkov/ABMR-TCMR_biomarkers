library("sva")
library("mygene")

# Helpful functions ------
cleanStudyData <- function(esetList, sampleMetadata) {
    ematList = foreach(studyName = names(esetList)) %do% {
        keepIdx = colnames(esetList[[studyName]]) %in% sampleMetadata[sampleMetadata$study == studyName, 'sample']
        exprs(esetList[[studyName]])[,keepIdx]
    }
    names(ematList) = names(esetList)
    return(ematList)
}

geneId2Symbol = function(emat) {
    geneSymbols = mygene::getGenes(rownames(emat), fields = "symbol", return.as = "DataFrame")
    rownames(emat) = geneSymbols$symbol
    
    cat("Unmapped", sum(is.na(geneSymbols$symbol)), "gene(s) with Id(s):", geneSymbols[is.na(geneSymbols$symbol), "query"], "\n")
    emat = emat[!is.na(rownames(emat)),]
    
    duplSymbols = rownames(emat)[duplicated(rownames(emat))]
    if (length(duplSymbols) > 0) {
        for (dupl in duplSymbols) {
            exprsTmp = emat[rownames(emat) %in% dupl,, drop = F]
            emat = emat[!(duplicated(rownames(emat)) & rownames(emat) == dupl),]
            emat[dupl,] = rowMeans(t(exprsTmp), na.rm = TRUE)
        }
    }
    return(emat)
}


# Data load, cleaning, combining -----
# Loading the metadata
# metadata = readxl::read_xlsx("../Metadata/Biopsy_Tx.xlsx", sheet = 1)

metadata = read.csv('../Metadata/kidney_GSM_metadata.csv', header = TRUE, stringsAsFactors = FALSE)
metadata = metadata[metadata$transplant == '1',]
metadata = metadata[,c(1,2,6)]
detect = str_detect(unname(unlist(metadata[,'subclass'])), "ABMR|TCMR|Mixed|BL|AR[+]CAN|AR|BL[+]CAN|STA|Normal")
metadata = metadata[detect,]

metadata[metadata$subclass %in% c("ABMR","TCMR", "Mixed","AR+CAN", "BL", "BL+CAN", "AR"), "class"] = "AR"
metadata[metadata$subclass == "STA", "class"] = "STA"
metadata[metadata$subclass == "Normal", "class"] = "Normal"
rownames(metadata) = metadata$gsm
colnames(metadata) = c("study", "sample", "subclass", "class")


# Loading all datasets
esetList = sapply(unique(metadata$study), function(x) readRDS(paste0("../Data_in/", x, ".rds")))

# Extracting expresision matrices
ematList = cleanStudyData(esetList, metadata)

# Reduce datasets to common genes
geneIds = Reduce(intersect, lapply(ematList, function(x) rownames(x)))
ematList1 = foreach(studyName = names(ematList)) %do% {ematNow = ematList[[studyName]][geneIds,]}

# Z-scale each dataset and combine them all into one matrix
ematMerged = do.call(cbind, ematList1)
emat.tmp = list()
for (batch in unique(metadata$batch)) {
    message("Processing batch ", batch, "...\n")
    samples = metadata[metadata$batch == batch, "sample"]
    tmp = ematMerged[, samples]
    tmp = (tmp - mean(tmp)) / sd(tmp)
    # tmp = normalizeQuantiles(tmp)
    emat.tmp[[batch]] = tmp
}
ematMerged = do.call(cbind, emat.tmp)

# Mapping entrez Ids to gene symbols
ematMerged = geneId2Symbol(ematMerged)

# Check match of column names with sample names from the metadata
all(colnames(ematMerged) == rownames(metadata))
ematMerged = ematMerged[,metadata$sample]

# Save non-batch corrected data
saveRDS(ematMerged, "../Data_out/data_merged_nonNorm.rds")
# ematMerged = readRDS("../Data_out/data_merged_nonNorm.rds")

# Batch correction -----
# Batch correction with preservation of bilogical differences in case of case-control unbalanced datasets
mod = model.matrix(~as.factor(class), data = metadata)
ematMerged.Combat = ComBat(ematMerged, 
                           batch = as.factor(metadata$batch), 
                           mod = mod,
                           # ref.batch = "GSE47683",
                           par.prior = TRUE, 
                           prior.plots = FALSE,
                           BPPARAM = MulticoreParam(3)) 

saveRDS(ematMerged.Combat, "../Data_out/data_merged_Norm.rds")